load('../System_Calibration/cmag_05-12-18.mat')

xv = linspace(-0.0525, 0.0525, 5);
yv = linspace(-0.04, 0.04, 5);
zv = linspace(-0.0475, 0.0875, 5);

[xd, yd, zd] = meshgrid(xv, yv, zv);

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

C = get_divfree_rbf_coefficients( nodes, values, eps);

%% Evaluation part
xv = linspace(-0.0525, 0.0525, 20);
yv = linspace(-0.04, 0.04, 20);
zv = linspace(-0.0475, 0.0875, 20);

[xd, yd, zd] = meshgrid(xv, yv, zv);

xd = reshape(xd, [], 1);
yd = reshape(yd, [], 1);
zd = reshape(zd, [], 1);

real = zeros(length(xd), 3);
interp = zeros(length(xd), 3);

for i = 1:length(xd)
    BG = cmag.FieldAndGradient([xd(i); yd(i); zd(i)], currents);
    real(i,:) = BG(1:3);
    interp(i,:) = evaluate_rbf([xd(i), yd(i), zd(i)], nodes, eps, C);
end

fprintf('average error: %f mT\n', 1000 * mean(sqrt(sum((real - interp).^2, 2))));
fprintf('median error: %f mT\n', 1000 * median(sqrt(sum((real - interp).^2, 2))));

figure;
hold on;
quiver3(xd, yd, zd, real(:,1), real(:,2), real(:,3));
quiver3(xd, yd, zd, interp(:,1), interp(:,2), interp(:,3));
hold off;
