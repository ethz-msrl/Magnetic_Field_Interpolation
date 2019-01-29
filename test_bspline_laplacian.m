%% Getting values
load('../System_Calibration/cmag_05-12-18.mat')

% discretization of input grid
nv = 5;
dim = 4;

xv = linspace(-0.0525, 0.0525, nv);
xstep = (xv(end) - xv(1)) / (length(xv) - 1);
yv = linspace(-0.04, 0.04, nv);
ystep = (yv(end) - yv(1)) / (length(yv) - 1);
zv = linspace(-0.0475, 0.0875, nv);
zstep = (zv(end) - zv(1)) / (length(zv) - 1);

[xd, yd, zd] = ndgrid(xv, yv, zv);

xd = xd(:);
yd = yd(:);
zd = zd(:);

nodes = [xd, yd, zd];

currents = [-1.963839324961165,-8.480666166183163,-5.201676928926839,-7.533621303296689,-6.321844234351666,-5.200949486701944,-1.654658618312610,-9.006911393485158, 1]';

values = zeros(length(xd), 3);
for i = 1:length(xd)
    BG = cmag.FieldAndGradient(nodes(i,:), currents);
    values(i, :) = BG(1:3);
end

%% Evaluation part

% normalize the values to the input grid
xve = linspace(-0.0525, 0.0525, 10);
yve = linspace(-0.04, 0.04, 10);
zve = linspace(-0.0475, 0.0875, 10);

[xde, yde, zde] = ndgrid(xve, yve, zve);

xde = reshape(xde, [], 1);
yde = reshape(yde, [], 1);
zde = reshape(zde, [], 1);

real = zeros(length(xde), 3);
for i = 1:length(xde)
    BG = cmag.FieldAndGradient([xde(i); yde(i); zde(i)], currents);
    real(i,:) = BG(1:3);
end

%% Spline Fit
nk = 5;

% this gets the appropriate knot sequence that satisfies the
% Schoenberg-Whitney condition
k_n = aptknt(xv, dim);
k_m = aptknt(yv, dim);
k_p = aptknt(zv, dim);

N = spcol(k_n, dim, xv);
M = spcol(k_m, dim, yv);
P = spcol(k_p, dim, zv);

Zi = kron(kron(N,M), P);
% we align the coefficients as [Cijkx Cijky; Cijkz]^T
Z = blkdiag(Zi, Zi, Zi);

% getting the first derivatives of the bases
temp = spcol(k_n, dim, brk2knt(xv, 2));
Ndot = temp(2:2:end, :);

temp = spcol(k_m, dim, brk2knt(yv, 2));
Mdot = temp(2:2:end, :);

temp = spcol(k_p, dim, brk2knt(zv, 2));
Pdot = temp(2:2:end, :);

% Zx = kron(kron(Ndot, M), P);
% Zy = kron(kron(N, Mdot), P);
% Zz = kron(kron(N, M), Pdot);

Zx = kron(P, kron(M, Ndot));
Zy = kron(P, kron(Mdot, N));
Zz = kron(Pdot, kron(M, N));

Q = [Zx, Zy, Zz; zeros(size(Zx)), -Zz, Zy; Zz, zeros(size(Zy)), -Zx; -Zy, Zx, zeros(size(Zz))];
%Q = [Zx, Zy, Zz];
D = reshape(values, size(values,1)*3, []);
C_ = Zi \ values;

C = quadprog(Z'*Z, -D'*Z, [], [], Q, zeros(size(Q,1),1));
%C = quadprog(Z'*Z, -D'*Z);

%% Check divergence
temp = reshape(C, 5, 5, 5, 3);
pp = fn2fm(spmak({k_n, k_m, k_p}, permute(temp, [4, 1, 2, 3])), 'pp');

dBdx = fnder(pp, [1,0,0]);
dBdy = fnder(pp, [0,1,0]);
dBdz = fnder(pp, [0,0,1]);

tempx = fnval(dBdx, nodes');
tempy = fnval(dBdy, nodes');
tempz = fnval(dBdz, nodes');

div = tempx(1,:) + tempy(2,:) + tempz(3,:);

fprintf('mean divergence: %s mT\n', 1000 * mean(div));

%% Check curl
cx = tempz(2,:) - tempy(3,:);
cy = tempx(3,:) - tempz(1,:);
cz = tempy(1,:) - tempx(2,:);

%% Compare

N = spcol(k_n, dim, xve);
M = spcol(k_m, dim, yve);
P = spcol(k_p, dim, zve);

temp = spcol(k_n, dim, brk2knt(xve, 2));
Ndot = temp(2:2:end, :);

temp = spcol(k_m, dim, brk2knt(yve, 2));
Mdot = temp(2:2:end, :);

temp = spcol(k_p, dim, brk2knt(zve, 2));
Pdot = temp(2:2:end, :);

Zx = kron(P, kron(M, Ndot));
Zy = kron(P, kron(Mdot, N));
Zz = kron(Pdot, kron(M, N));

Q = [Zx, Zy, Zz];
fprintf('mean divergence: %s mT/m\n', 1000*mean(Q * C));
fprintf('max divergence: %s mT/m\n', 1000*max(Q * C));

Cx = [zeros(size(Zx)), -Zz, Zy];
Cy = [Zz, zeros(size(Zy)), -Zx];
Cz = [-Zy, Zx, zeros(size(Zz))];
fprintf('mean curl magnitude %s mT/m\n', 1000*mean(sqrt(sum([Cx * C, Cy * C, Cz * C].^2, 2))));
fprintf('max curl magnitude %s mT/m\n', 1000*max(sqrt(sum([Cx * C, Cy * C, Cz * C].^2, 2))));

Z = kron(kron(N,M), P);
interp = Z * reshape(C, size(values, 1), 3);

fprintf('average error: %f mT\n', 1000 * mean(sqrt(sum((real - interp).^2, 2))));
fprintf('median error: %f mT\n', 1000 * median(sqrt(sum((real - interp).^2, 2))));

figure;
hold on;
quiver3(xd, yd, zd, values(:,1), values(:,2), values(:,3));
quiver3(xde, yde, zde, real(:,1), real(:,2), real(:,3));
quiver3(xde, yde, zde, interp(:,1), interp(:,2), interp(:,3));
hold off;