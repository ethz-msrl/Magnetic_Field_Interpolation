%% Getting values
load('../../System_Calibration/cmag_05-12-18.mat')

% discretization of input grid
nv = 5;

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

D = 5; % degree of b-spline
nk = nv;

% this gets the appropriate knot sequence that satisfies the
% Schoenberg-Whitney condition
k_n = aptknt(linspace(xv(1), xv(end), nk), D);
k_m = aptknt(linspace(yv(1), yv(end), nk), D);
k_p = aptknt(linspace(zv(1), zv(end), nk), D);

N = spcol(k_n, D, xv);
M = spcol(k_m, D, yv);
P = spcol(k_p, D, zv);

Z = kron(kron(N,M), P);
C = Z \ values;

%% Compare

N = spcol(k_n, D, xve);
M = spcol(k_m, D, yve);
P = spcol(k_p, D, zve);

Z = kron(kron(N,M), P);
interp = Z * C;

fprintf('average error: %f mT\n', 1000 * mean(sqrt(sum((real - interp).^2, 2))));
fprintf('median error: %f mT\n', 1000 * median(sqrt(sum((real - interp).^2, 2))));

figure;
hold on;
quiver3(xd, yd, zd, values(:,1), values(:,2), values(:,3));
quiver3(xde, yde, zde, real(:,1), real(:,2), real(:,3));
quiver3(xde, yde, zde, interp(:,1), interp(:,2), interp(:,3));
hold off;

%% Kernel size

% we iterate the number of knots and see the average fit error
nkv = 5:16;
avg_errors = zeros(size(nkv));

for i = 1:length(nkv)

    nk = nkv(i);
    
%      k_n = optknt(linspace(xv(1), xv(end), nk), D, 100);
%      k_m = optknt(linspace(yv(1), yv(end), nk), D, 100);
%      k_p = optknt(linspace(zv(1), zv(end), nk), D, 100);

    k_n = aptknt(linspace(xv(1), xv(end), nk), D);
    k_m = aptknt(linspace(yv(1), yv(end), nk), D);
    k_p = aptknt(linspace(zv(1), zv(end), nk), D);
%     k_n = augknt(linspace(xv(1), xv(end), nk), D);
%     k_m = augknt(linspace(yv(1), yv(end), nk), D);
%     k_p = augknt(linspace(zv(1), zv(end), nk), D);

    N = spcol(k_n, D, xv);
    M = spcol(k_m, D, yv);
    P = spcol(k_p, D, zv);

    Z = kron(kron(N,M), P);
    C = Z \ values;

    N = spcol(k_n, D, xve);
    M = spcol(k_m, D, yve);
    P = spcol(k_p, D, zve);

    Z = kron(kron(N,M), P);
    interp = Z * C;

    avg_errors(i) = 1000 * mean(sqrt(sum((real - interp).^2, 2)));

end

bar(nkv, avg_errors);

%% Grid size
nvv = 3:8;
avg_errors = zeros(size(nvv));
nk = 4;

for i = 1:length(nvv)
    nv = nvv(i);
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
    for j = 1:length(xd)
        BG = cmag.FieldAndGradient(nodes(j,:), currents);
        values(j, :) = BG(1:3);
    end
    
    k_n = aptknt(linspace(xv(1), xv(end), nk), D);
    k_m = aptknt(linspace(yv(1), yv(end), nk), D);
    k_p = aptknt(linspace(zv(1), zv(end), nk), D);

    N = spcol(k_n, D, xv);
    M = spcol(k_m, D, yv);
    P = spcol(k_p, D, zv);

    Z = kron(kron(N,M), P);
    C = Z \ values;
    
    N = spcol(k_n, D, xve);
    M = spcol(k_m, D, yve);
    P = spcol(k_p, D, zve);

    Z = kron(kron(N,M), P);
    interp = Z * C;

    avg_errors(i) = 1000 * mean(sqrt(sum((real - interp).^2, 2)));
end

figure;
bar(nvv, avg_errors);
