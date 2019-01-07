%% Interpolation part
[M, B_fun] = get_tricubic_div_free_matrix();

syms x y z

SCALE = 30e-3;

zv = SCALE * [0, 1, 0, 1, 0, 1, 0, 1];
yv = SCALE * [0, 0, 1, 1, 0, 0, 1, 1];
xv = SCALE * [0, 0, 0, 0, 1, 1, 1, 1];

currents = [6.7573    4.3353   -7.3174   -1.9490    3.2694    3.6722   -4.4116   -3.6951]';

fields = zeros(8, 3);
grads = zeros(8,9);

Bx_x = zeros(1,8);
By_y = zeros(1,8);
Bx_y = zeros(1,8);
Bx_z = zeros(1,8);
By_z = zeros(1,8);

eps = SCALE;

for i = 1:8
    pos = [xv(i), yv(i), zv(i)]; %+ [SCALE,0,0];
    BG = magSystem.FieldAndGradient(pos, currents);
    fields(i,:) = BG(1:3);
    
    pos_px = [xv(i) + eps, yv(i), zv(i)];
    pos_nx = [xv(i) - eps, yv(i), zv(i)];
    BG_px = magSystem.FieldAndGradient(pos_px, currents);
    BG_nx = magSystem.FieldAndGradient(pos_nx, currents);
    
    pos_py = [xv(i), yv(i) + eps, zv(i)];
    pos_ny = [xv(i), yv(i) - eps, zv(i)];
    BG_py = magSystem.FieldAndGradient(pos_py, currents);
    BG_ny = magSystem.FieldAndGradient(pos_ny, currents);
    
    pos_pz = [xv(i), yv(i), zv(i) + eps];
    pos_nz = [xv(i), yv(i), zv(i) - eps];
    BG_pz = magSystem.FieldAndGradient(pos_pz, currents);
    BG_nz = magSystem.FieldAndGradient(pos_nz, currents);
    
    pos_py_pz = [xv(i), yv(i) + eps, zv(i) + eps];
    pos_py_nz = [xv(i), yv(i) + eps, zv(i) - eps];
    pos_ny_pz = [xv(i), yv(i) - eps, zv(i) + eps];
    pos_ny_nz = [xv(i), yv(i) - eps, zv(i) - eps];
    
    BG_py_pz = magSystem.FieldAndGradient(pos_py_pz, currents);
    BG_py_nz = magSystem.FieldAndGradient(pos_py_nz, currents);
    BG_ny_pz = magSystem.FieldAndGradient(pos_ny_pz, currents);
    BG_ny_nz = magSystem.FieldAndGradient(pos_ny_nz, currents);
    
    % dBx/dx
    Bx_x(:) = (BG_px(1) - BG_nx(1))/(2);
    grads(i,1) = Bx_x;
    
    % dBy/dy
    By_y(:) = (BG_py(2) - BG_ny(2))/(2);
    grads(i,5) = By_y;
    
    % dBx/dy
    Bx_y(:) = (BG_py(1) - BG_ny(1))/(2);
    grads(i,2) = Bx_y;
    grads(i,4) = Bx_y; 
    
    % dBx/dz
    Bx_z(:) = (BG_pz(1) - BG_nz(1))/2;
    grads(i,3) = Bx_z;
    grads(i,7) = Bx_z;
    
    % dBy/dz
    By_z(:) = (BG_pz(2) - BG_nz(2))/2;
    grads(i,6) = By_z;
    grads(i,8) = By_z;
    
    grads(i, 9) = Bx_x + By_y;
end

% Calculate derivatives using forward differences
% eps = 1e-8;
% Bx_x = zeros(1,8);
% By_y = zeros(1,8);
% Bx_y = zeros(1,8);
% Bx_z = zeros(1,8);
% By_z = zeros(1,8);
% Bx_yz = zeros(1,8);
% 
% for i = 1:8
%     pos_px = [xv(i) + eps, yv(i), zv(i)];
%     pos_nx = [xv(i) - eps, yv(i), zv(i)];
%     BG_px = magSystem.FieldAndGradient(pos_px, currents);
%     BG_nx = magSystem.FieldAndGradient(pos_nx, currents);
%     
%     pos_py = [xv(i), yv(i) + eps, zv(i)];
%     pos_ny = [xv(i), yv(i) - eps, zv(i)];
%     BG_py = magSystem.FieldAndGradient(pos_py, currents);
%     BG_ny = magSystem.FieldAndGradient(pos_ny, currents);
%     
%     pos_pz = [xv(i), yv(i), zv(i) + eps];
%     pos_nz = [xv(i), yv(i), zv(i) - eps];
%     BG_pz = magSystem.FieldAndGradient(pos_pz, currents);
%     BG_nz = magSystem.FieldAndGradient(pos_nz, currents);
%     
%     pos_py_pz = [xv(i), yv(i) + eps, zv(i) + eps];
%     pos_py_nz = [xv(i), yv(i) + eps, zv(i) - eps];
%     pos_ny_pz = [xv(i), yv(i) - eps, zv(i) + eps];
%     pos_ny_nz = [xv(i), yv(i) - eps, zv(i) - eps];
%     
%     BG_py_pz = magSystem.FieldAndGradient(pos_py_pz, currents);
%     BG_py_nz = magSystem.FieldAndGradient(pos_py_nz, currents);
%     BG_ny_pz = magSystem.FieldAndGradient(pos_ny_pz, currents);
%     BG_ny_nz = magSystem.FieldAndGradient(pos_ny_nz, currents);
%     
%     % dBx/dx
%     Bx_x(:) = (BG_px(1) - BG_nx(1))/(2);
%     
%     % dBy/dy
%     By_y(:) = (BG_py(2) - BG_ny(2))/(2);
%     
%     % dBx/dy
%     Bx_y(:) = (BG_py(1) - BG_ny(1))/(2);
%     
%     % dBx/dz
%     Bx_z(:) = (BG_pz(1) - BG_nz(1))/2;
%     
%     % dBy/dz
%     By_z(:) = (BG_pz(2) - BG_nz(2))/2;
%     
%     % d2Bx/dydz
%     Bx_yz(:,i) = (BG_py_pz(1) - BG_ny_pz(1) - BG_py_nz(1) + BG_ny_nz(1)) / (4);
% end

D = zeros(48,1);
D(1:24) = [
    reshape(fields', 8*3,1);
    ];
%     Bx_x';
%     Bx_y';
%     Bx_z';
%     By_y';
%     By_z'];

a_sol = M \ D;

% test in plane
N = 8;

xvt = linspace(0,1,N);
yvt = linspace(0,1,N);

[xvtg, yvtg] = meshgrid(xvt, yvt);

fields_interp = zeros(N, N, 3);
fields_valid = zeros(N, N, 3);

for i = 1:N
    for j = 1:N
        pos = SCALE * [xvtg(i,j), yvtg(i,j), 0];
        fields_interp(i, j, :) = B_fun(xvtg(i, j), yvtg(i,j), 0, a_sol);
        BG = magSystem.FieldAndGradient(pos, currents);
        fields_valid(i,j,:) = BG(1:3);
    end
end

q1 = quiver(xvtg, yvtg, squeeze(fields_interp(:,:,1)), ...
    squeeze(fields_interp(:,:,2)), 0);
hold on;
q2 = quiver(xvtg, yvtg, squeeze(fields_valid(:,:,1)), squeeze(fields_valid(:,:,2)), 0);
hold off;

scale = 5;
qU1 = get(q1, 'UData');
qV1 = get(q1, 'VData');
set(q1, 'UData', scale*qU1, 'VData', scale*qV1);
qU2 = get(q2, 'UData');
qV2 = get(q2, 'VData');
set(q2, 'UData', scale*qU2, 'VData', scale*qV2);

% for 3D Data

% xvt = linspace(0,1,N);
% yvt = linspace(0,1,N);
% zvt = linspace(0,1,N);
% 
% [xvtg, yvtg, zvtg] = meshgrid(xvt, yvt, zvt);
% 
% fields_interp = zeros(N, N, N, 3);
% fields_valid = zeros(N, N, N, 3);
% 
% for i = 1:N
%     for j = 1:N
%         for k =1:N
%             pos = SCALE* [xvt(i), yvt(j), zvt(k)];
%             fields_interp(i, j, k, :) = B_fun(xvt(i)*SCALE, yvt(j)*SCALE, zvt(k));
%             BG = magSystem.FieldAndGradient(pos, currents);
%             fields_valid(i,j,k,:) = BG(1:3);
%         end
%     end
% end
% 
% q1 = quiver3(xvtg, yvtg, zvtg, fields_interp(:,:,:,1), fields_interp(:,:,:,2),...
%     fields_interp(:,:,:,3));
% hold on;
% q2 = quiver3(xv/SCALE, yv/SCALE, zv/SCALE, fields(:,1)', fields(:,2)', ...
%     fields(:,3)');
% q3 = quiver3(xvtg, yvtg, zvtg, fields_valid(:,:,:,1), fields_valid(:,:,:,2),...
%     fields_valid(:,:,:,3));
% hold off;
        