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

syms x y z

%% Evaluation part

% normalize the values to the input grid
xve = linspace(-0.0525, 0.0525, 10);
yve = linspace(-0.04, 0.04, 10);
zve = linspace(-0.0475, 0.0875, 10);

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

for i = 1:length(xd)
    ix = floor(xd(i)) + 1;
    iy = floor(yd(i)) + 1;
    iz = floor(zd(i)) + 1;
    
    x = xd(i) + 1 - ix;
    y = yd(i) + 1 - iy;
    z = zd(i) + 1 - iz;
    
    % handling points that are on the extremities of the data,
    % we don't use find the next cell that contains but the previous
    if ix == (xven(end) + 1)
        ix = ix - 1;
        x = 1;
    end
    
    if iy == (yven(end) + 1)
        iy = iy - 1;
        y = 1;
    end
    
    if iz == (zven(end) + 1)
        iz = iz - 1;
        z = 1;
    end
    
    fields = [
        squeeze(vg(ix, iy, iz, :));
        squeeze(vg(ix, iy, iz+1, :));
        squeeze(vg(ix, iy+1, iz, :));
        squeeze(vg(ix, iy+1, iz+1, :));
        squeeze(vg(ix+1, iy, iz, :));
        squeeze(vg(ix+1, iy, iz+1, :));
        squeeze(vg(ix+1, iy+1, iz, :));
        squeeze(vg(ix+1, iy+1, iz+1, :));
        ];
    
    interp(i, :) = evaluate_tricubic_div_free(fields, x, y, z, M, B_fun);
end

fprintf('average error: %f mT\n', 1000 * mean(sqrt(sum((real - interp).^2, 2))));

figure;
hold on;
quiver3(xd, yd, zd, real(:,1), real(:,2), real(:,3));
quiver3(xd, yd, zd, interp(:,1), interp(:,2), interp(:,3));
hold off;

