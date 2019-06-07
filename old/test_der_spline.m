% A small test script to check that the you can take the derivative 
% of 2D B-Splines by using the collocation derivatives

dim = 4;
xv = linspace(-1, 1, 5);
yv = linspace(-1.1, 1.1, 5);

k_n = aptknt(xv, dim);
k_m = aptknt(yv, dim);

N = spcol(k_n, dim, xv);
temp = spcol(k_n, dim, brk2knt(xv, 2));
Ndot = temp(2:2:end, :);
M = spcol(k_m, dim, yv);
temp = spcol(k_m, dim, brk2knt(yv, 2));
Mdot = temp(2:2:end, :);

[X, Y] = ndgrid(xv, yv);
values = exp(-2.*(X.^2 + Y.^2));

C = kron(N,M) \ reshape(values, length(xv)*length(yv), []);

temp = reshape(C, 5, 5);
f = spmak({k_n, k_m}, temp);

dfdx = fnder(f, [1,0]);

xve = linspace(xv(1), xv(end), 101);
yve = linspace(xv(1), yv(end), 101);

v_f = fnval(f, {xve, yve});
v_dfdx = fnval(dfdx, {xve, yve});

subplot(3,1,1);
surf(v_f);
xlabel('x');

subplot(3,1,2);
surf(v_dfdx);
xlabel('x');

Me = spcol(k_m, dim, xve);
Ne = spcol(k_n, dim, yve);
temp = spcol(k_n, dim, brk2knt(xve, 2));
Ndote = temp(2:2:end, :);
temp = spcol(k_m, dim, brk2knt(yve, 2));
Mdote = temp(2:2:end, :);
v_dfdx_ = reshape(kron(Ndote, Me) * C, 101, 101)';

% the first plot is the function value
% the second is the x-derivative calculated by taking the derivative
% of the b-form function
% the third is the x-derivative by calculating the collocation points and
% using the coefficients
% the third plot should match the second one

subplot(3,1,3);
surf(v_dfdx_);
xlabel('x');

%% now for a 2D incompressible vector field
xve = linspace(xv(1), xv(end), 16);
yve = linspace(yv(1), yv(end), 16);
[Xe, Ye] = ndgrid(xve, yve);

N = spcol(k_n, dim, xv);
temp = spcol(k_n, dim, brk2knt(xv, 2));
Ndot = temp(2:2:end, :);
M = spcol(k_m, dim, yv);
temp = spcol(k_m, dim, brk2knt(yv, 2));
Mdot = temp(2:2:end, :);
Ne = spcol(k_n, dim, xve);
Me = spcol(k_m, dim, yve);
temp = spcol(k_n, dim, brk2knt(xve, 2));
Ndote = temp(2:2:end, :);
temp = spcol(k_m, dim, brk2knt(yve, 2));
Mdote = temp(2:2:end, :);

values = cat(3, -Y, X);
D = reshape(values, length(xv)*length(yv)*2, []);
Zi = kron(N, M);
Z = blkdiag(Zi, Zi);

%Zx = kron(Ndot, M);
%Zy = kron(N, Mdot);
% WTF do I need to inverse the order of the kron product to get the 
% right partial derivative?
Zx = kron(M, Ndot);
Zy = kron(Ndot, M);

Q = [Zx, Zy];

%C = quadprog(Z'*Z, -D'*Z, [], [], Q, zeros(size(Q,1),1));
C = quadprog(Z'*Z, -D'*Z);
C = reshape(C, length(xv)*length(yv), 2);
% C = Z \ reshape(values, length(xv)*length(yv), 2);

temp = reshape(C, 5, 5, 2);
temp = permute(temp, [3, 1, 2]);
f = spmak({k_n, k_m}, temp);
v_f = fnval(f, {xve, yve});

figure;

quiver(Xe, Ye, squeeze(v_f(1,:,:)), squeeze(v_f(2,:,:)));
hold on;
quiver(X, Y, values(:,:,1), values(:,:,2));

df_dx = fnder(f, [1, 0]);
df_dy = fnder(f, [0, 1]);

temp1 = fnval(df_dx, {xve, yve});
temp2 = fnval(df_dy, {xve, yve});

div = squeeze(temp1(1,:,:) + temp2(2,:,:));
%surf(Xe, Ye, div)

%temp = kron(Ndote, Me) * C;
temp = kron(Me, Ndote) * C;
vx_dx = reshape(temp(:,1), 16, 16);
%temp = kron(Ne, Mdote) * C;
temp = kron(Mdote, Ne) * C;
vy_dy = reshape(temp(:,2), 16, 16);

div_ = vx_dx + vy_dy;

fprintf('avg divergence error: %f\n', mean(mean(abs(div - div_))));

