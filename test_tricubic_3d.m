xv = linspace(-0.0525, 0.0525, 5);
xstep = (xv(end) - xv(1)) / (length(xv) - 1);
yv = linspace(-0.04, 0.04, 5);
ystep = (yv(end) - yv(1)) / (length(yv) - 1);
zv = linspace(-0.0475, 0.0875, 5);
zstep = (zv(end) - zv(1)) / (length(zv) - 1);

Nx = length(xv);
Ny = length(yv);
Nz = length(zv);

[xd, yd, zd] = ndgrid(xv, yv, zv);

xd = xd(:);
yd = yd(:);
zd = zd(:);

nodes = [xd, yd, zd];

currents = [-1.963839324961165,-8.480666166183163,-5.201676928926839,-7.533621303296689,-6.321844234351666,-5.200949486701944,-1.654658618312610,-9.006911393485158, 1]';

values = zeros(length(xd), 3);
grad = zeros(length(xd), 5);
for i = 1:length(xd)
    BG = cmag.FieldAndGradient(nodes(i,:), currents);
    values(i, :) = BG(1:3);
    grad(i,:) = BG(4:end);
end

grad = reshape(grad, [Nx, Ny, Nz, 5]);

[M, B_fun, G_fun] = get_tricubic_3d_matrix();

syms x y z

%% Evaluation part

% normalize the values to the input grid
xve = linspace(-0.0525, 0.0524, 10);
yve = linspace(-0.04, 0.039, 10);
zve = linspace(-0.0475, 0.0874, 10);

[xd, yd, zd] = ndgrid(xve, yve, zve);

xd = reshape(xd, [], 1);
yd = reshape(yd, [], 1);
zd = reshape(zd, [], 1);

real = zeros(length(xd), 3);

for i = 1:length(xd)
    BG = cmag.FieldAndGradient([xd(i); yd(i); zd(i)], currents);
    real(i,:) = BG(1:3);
end

%% Calculation Part

xven = (xve - xv(1)) / xstep;
yven = (yve - yv(1)) / ystep;
zven = (zve - zv(1)) / zstep;

[xd, yd, zd] = ndgrid(xven, yven, zven);

xd = reshape(xd, [], 1);
yd = reshape(yd, [], 1);
zd = reshape(zd, [], 1);
% 
interp = zeros(length(xd), 3);

vg = reshape(values, [Nx, Ny, Nz, 3]);

[dvg_dy, dvg_dx, vg_dz] = gradient(vg);
[dvg_dxy, ~, dvg_dxz] = gradient(dvg_dx);
[~,~,dvg_dyz] = gradient(dvg_dy);
[~,~,dvg_dxyz] = gradient(dvg_dxy);


for i = 1:length(xd)

    ix = floor(xd(i)) + 1;
    iy = floor(yd(i)) + 1;
    iz = floor(zd(i)) + 1;
    
    x = xd(i) + 1 - ix;
    y = yd(i) + 1 - iy;
    z = zd(i) + 1 - iz;
   
    D = squeeze([
    % f
        vg(ix, iy, iz, :);
        vg(ix, iy, iz+1, :);
        vg(ix, iy+1, iz, :);
        vg(ix, iy+1, iz+1, :);
        vg(ix+1, iy, iz, :);
        vg(ix+1, iy, iz+1, :);
        vg(ix+1, iy+1, iz, :);
        vg(ix+1, iy+1, iz+1, :);
    % df/dx
        dvg_dx(ix, iy, iz, :);
        dvg_dx(ix, iy, iz+1, :);
        dvg_dx(ix, iy+1, iz, :);
        dvg_dx(ix, iy+1, iz+1, :);
        dvg_dx(ix+1, iy, iz, :);
        dvg_dx(ix+1, iy, iz+1, :);
        dvg_dx(ix+1, iy+1, iz, :);
        dvg_dx(ix+1, iy+1, iz+1, :);
    % df/dy
        dvg_dy(ix, iy, iz, :);
        dvg_dy(ix, iy, iz+1, :);
        dvg_dy(ix, iy+1, iz, :);
        dvg_dy(ix, iy+1, iz+1, :);
        dvg_dy(ix+1, iy, iz, :);
        dvg_dy(ix+1, iy, iz+1, :);
        dvg_dy(ix+1, iy+1, iz, :);
        dvg_dy(ix+1, iy+1, iz+1, :);
    % df/dz
        dvg_dz(ix, iy, iz, :);
        dvg_dz(ix, iy, iz+1, :);
        dvg_dz(ix, iy+1, iz, :);
        dvg_dz(ix, iy+1, iz+1, :);
        dvg_dz(ix+1, iy, iz, :);
        dvg_dz(ix+1, iy, iz+1, :);
        dvg_dz(ix+1, iy+1, iz, :);
        dvg_dz(ix+1, iy+1, iz+1, :);
    % d2f/dxdy
        dvg_dxy(ix, iy, iz, :);
        dvg_dxy(ix, iy, iz+1, :);
        dvg_dxy(ix, iy+1, iz, :);
        dvg_dxy(ix, iy+1, iz+1, :);
        dvg_dxy(ix+1, iy, iz, :);
        dvg_dxy(ix+1, iy, iz+1, :);
        dvg_dxy(ix+1, iy+1, iz, :);
        dvg_dxy(ix+1, iy+1, iz+1, :);
    % d2f/dxdz
        dvg_dxz(ix, iy, iz, :);
        dvg_dxz(ix, iy, iz+1, :);
        dvg_dxz(ix, iy+1, iz, :);
        dvg_dxz(ix, iy+1, iz+1, :);
        dvg_dxz(ix+1, iy, iz, :);
        dvg_dxz(ix+1, iy, iz+1, :);
        dvg_dxz(ix+1, iy+1, iz, :);
        dvg_dxz(ix+1, iy+1, iz+1, :);
     % d2f/dydz
        dvg_dyz(ix, iy, iz, :);
        dvg_dyz(ix, iy, iz+1, :);
        dvg_dyz(ix, iy+1, iz, :);
        dvg_dyz(ix, iy+1, iz+1, :);
        dvg_dyz(ix+1, iy, iz, :);
        dvg_dyz(ix+1, iy, iz+1, :);
        dvg_dyz(ix+1, iy+1, iz, :);
        dvg_dyz(ix+1, iy+1, iz+1, :);
     % d3f/dxdydz
        dvg_dxyz(ix, iy, iz, :);
        dvg_dxyz(ix, iy, iz+1, :);
        dvg_dxyz(ix, iy+1, iz, :);
        dvg_dxyz(ix, iy+1, iz+1, :);
        dvg_dxyz(ix+1, iy, iz, :);
        dvg_dxyz(ix+1, iy, iz+1, :);
        dvg_dxyz(ix+1, iy+1, iz, :);
        dvg_dxyz(ix+1, iy+1, iz+1, :);
     ]);
 
    a_sol = M \ D;

     interp(i, :) = [  B_fun(x, y, z, a_sol(:,1)'), ...
            B_fun(x, y, z, a_sol(:,2)'), ...
            B_fun(x, y, z, a_sol(:,3)')];
end

fprintf('average error: %f mT\n', 1000 * mean(sqrt(sum((real - interp).^2, 2))));
fprintf('median error: %f mT\n', 1000 * median(sqrt(sum((real - interp).^2, 2))));

figure;
hold on;
quiver3(xd, yd, zd, real(:,1), real(:,2), real(:,3));
quiver3(xd, yd, zd, interp(:,1), interp(:,2), interp(:,3));
hold off;

