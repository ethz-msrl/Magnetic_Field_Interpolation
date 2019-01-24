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

currents = [-1.963839324961165,-8.480666166183163,-5.201676928926839,-7.533621303296689,-6.321844234351666,-5.200949486701944,-1.654658618312610,-9.006911393485158, 1]';

values = zeros(3, length(xv), length(yv), length(zv));
% values = zeros(length(xd), 3);
% for i = 1:length(xd)
%     BG = cmag.FieldAndGradient(nodes(i,:), currents);
%     values(i, :) = BG(1:3);
% end

for i = 1:length(xv)
    for j = 1:length(yv)
        for k = 1:length(zv)
            BG = cmag.FieldAndGradient([xv(i), yv(j), zv(k)], currents);
            values(:, i, j, k) = BG(1:3);
        end
    end
end

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

%% Compare
pp = csapi({xv, yv, zv}, values);
interp = fnval(pp, [xd, yd, zd]')';

fprintf('average error: %f mT\n', 1000 * mean(sqrt(sum((real - interp).^2, 2))));
fprintf('median error: %f mT\n', 1000 * median(sqrt(sum((real - interp).^2, 2))));

figure;
hold on;
quiver3(xd, yd, zd, real(:,1), real(:,2), real(:,3));
quiver3(xd, yd, zd, interp(:,1), interp(:,2), interp(:,3));
hold off;
