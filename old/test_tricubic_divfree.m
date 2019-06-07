xv = linspace(-0.0525, 0.0525, 5);
xstep = (xv(end) - xv(1)) / (length(xv) - 1);
yv = linspace(-0.04, 0.04, 5);
ystep = (yv(end) - yv(1)) / (length(yv) - 1);
zv = linspace(-0.0475, 0.0875, 5);
zstep = (zv(end) - zv(1)) / (length(zv) - 1);

[xd, yd, zd] = ndgrid(xv, yv, zv);

xd = xd(:);
yd = yd(:);
zd = zd(:);

nodes = [xd, yd, zd];

eps = 21;

currents = [-1.963839324961165,-8.480666166183163,-5.201676928926839,-7.533621303296689,-6.321844234351666,-5.200949486701944,-1.654658618312610,-9.006911393485158, 1]';

values = zeros(length(xd), 3);
for i = 1:length(xd)
    BG = cmag.FieldAndGradient(nodes(i,:), currents);
    values(i, :) = BG(1:3);
end

[M, B_fun] = get_tricubic_div_free_matrix();


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

xven = (xve - xv(1)) / xstep;
yven = (yve - yv(1)) / ystep;
zven = (zve - zv(1)) / zstep;

[xd, yd, zd] = ndgrid(xven, yven, zven);

xd = reshape(xd, [], 1);
yd = reshape(yd, [], 1);
zd = reshape(zd, [], 1);

interp = zeros(length(xd), 3);

Nx = length(xv);
Ny = length(yv);
Nz = length(zv);
vg = reshape(values, [Nx, Ny, Nz, 3]);
% we padd the array so we can do i-1 and i+2 on the borders
vg = padarray(vg, [2, 2, 2, 0], 'circular', 'both');

for i = 1:length(xd)
    % add 2 for the pad on left 
    ix = floor(xd(i)) + 3;
    iy = floor(yd(i)) + 3;
    iz = floor(zd(i)) + 3;
    
    x = xd(i) + 3 - ix;
    y = yd(i) + 3 - iy;
    z = zd(i) + 3 - iz;
    
    % dB/dx
    B_x = squeeze([
    0.5 * (vg(ix+1, iy, iz, :)       - vg(ix-1, iy, iz, :));
    0.5 * (vg(ix+1, iy, iz+1, :)     - vg(ix-1, iy, iz+1, :));
    0.5 * (vg(ix+1, iy+1, iz, :)     - vg(ix-1, iy+1, iz, :));
    0.5 * (vg(ix+1, iy+1, iz+1, :)   - vg(ix-1, iy+1, iz+1, :));
    0.5 * (vg(ix+2, iy, iz, :)       - vg(ix, iy, iz, :));
    0.5 * (vg(ix+2, iy, iz+1, :)     - vg(ix, iy, iz+1, :));
    0.5 * (vg(ix+2, iy+1, iz, :)     - vg(ix, iy+1, iz, :));
    0.5 * (vg(ix+2, iy+1, iz+1, :)   - vg(ix, iy+1, iz+1, :));
    ]);
    
    % dB/dy
    B_y = squeeze([0.5 * (vg(ix, iy+1, iz, :)       - vg(ix, iy-1, iz, :));
    0.5 * (vg(ix, iy+1, iz+1, :)     - vg(ix, iy-1, iz+1, :));
    0.5 * (vg(ix, iy+2, iz, :)       - vg(ix, iy, iz, :));
    0.5 * (vg(ix, iy+2, iz+1, :)     - vg(ix, iy, iz+1, :));
    0.5 * (vg(ix+1, iy+1, iz, :)     - vg(ix+1, iy-1, iz, :));
    0.5 * (vg(ix+1, iy+1, iz+1, :)   - vg(ix+1, iy-1, iz+1, :));
    0.5 * (vg(ix+1, iy+2, iz, :)     - vg(ix+1, iy, iz, :));
    0.5 * (vg(ix+1, iy+2, iz+1, :)   - vg(ix+1, iy, iz+1, :));
    ]);
    
    % dB/dz
    B_z = squeeze([0.5 * (vg(ix, iy, iz+1, :)       - vg(ix, iy, iz-1, :));
    0.5 * (vg(ix, iy, iz+2, :)       - vg(ix, iy, iz, :));
    0.5 * (vg(ix, iy+1, iz+1, :)     - vg(ix, iy+1, iz-1, :));
    0.5 * (vg(ix, iy+1, iz+2, :)     - vg(ix, iy+1, iz, :));
    0.5 * (vg(ix+1, iy, iz+1, :)     - vg(ix+1, iy, iz-1, :));
    0.5 * (vg(ix+1, iy, iz+2, :)     - vg(ix+1, iy, iz, :));
    0.5 * (vg(ix+1, iy+1, iz+1, :)   - vg(ix+1, iy+1, iz-1, :));
    0.5 * (vg(ix+1, iy+1, iz+2, :)   - vg(ix+1, iy+1, iz, :));
    ]);
    
    f = squeeze([
    % f
        vg(ix, iy, iz, :);
        vg(ix, iy, iz+1, :);
        vg(ix, iy+1, iz, :);
        vg(ix, iy+1, iz+1, :);
        vg(ix+1, iy, iz, :);
        vg(ix+1, iy, iz+1, :);
        vg(ix+1, iy+1, iz, :);
        vg(ix+1, iy+1, iz+1, :);
        ]);
    
    D = [reshape(f', [8*3, 1]); 
        reshape([B_x(:,1), B_y(:,1), B_z(:,1), B_x(:,2), B_y(:,2), B_z(:,2), B_x(:,3), B_y(:,3), B_z(:,3) ]', [8*9, 1])
        ];
    
    a_sol = M \ D;
    interp(i, :) = B_fun(x, y, z, a_sol);
end

fprintf('average error: %f mT\n', 1000 * mean(sqrt(sum((real - interp).^2, 2))));
fprintf('median error: %f mT\n', 1000 * median(sqrt(sum((real - interp).^2, 2))));

figure;
hold on;
quiver3(xd, yd, zd, real(:,1), real(:,2), real(:,3));
quiver3(xd, yd, zd, interp(:,1), interp(:,2), interp(:,3));
hold off;

