syms x y z

%% Interpolation part
SCALE = 30e-3;

zv = SCALE * [0, 1, 0, 1, 0, 1, 0, 1];
yv = SCALE * [0, 0, 1, 1, 0, 0, 1, 1];
xv = SCALE * [0, 0, 0, 0, 1, 1, 1, 1];

%currents = 16 * rand(8,1) - 8;
currents = [6.7573    4.3353   -7.3174   -1.9490    3.2694    3.6722   -4.4116   -3.6951]';

fields = zeros(8, 3);

for i = 1:8
    pos = [xv(i), yv(i), zv(i)]; %+ [SCALE,0,0];
    BG = cmag.FieldAndGradient(pos, currents);
    fields(i,:) = BG(1:3);
end

% Calculate derivatives using forward differences
eps = 1e-8;
Bx_y = zeros(1,8);
Bx_z = zeros(1,8);
By_z = zeros(1,8);
Bx_yz = zeros(1,8);

for i = 1:8
    pos_px = [xv(i) + eps, yv(i), zv(i)];
    pos_nx = [xv(i) - eps, yv(i), zv(i)];
    BG_px = cmag.FieldAndGradient(pos_px, currents);
    BG_nx = cmag.FieldAndGradient(pos_nx, currents);
    
    pos_py = [xv(i), yv(i) + eps, zv(i)];
    pos_ny = [xv(i), yv(i) - eps, zv(i)];
    BG_py = cmag.FieldAndGradient(pos_py, currents);
    BG_ny = cmag.FieldAndGradient(pos_ny, currents);
    
    pos_pz = [xv(i), yv(i), zv(i) + eps];
    pos_nz = [xv(i), yv(i), zv(i) - eps];
    BG_pz = cmag.FieldAndGradient(pos_pz, currents);
    BG_nz = cmag.FieldAndGradient(pos_nz, currents);
    
    pos_py_pz = [xv(i), yv(i) + eps, zv(i) + eps];
    pos_py_nz = [xv(i), yv(i) + eps, zv(i) - eps];
    pos_ny_pz = [xv(i), yv(i) - eps, zv(i) + eps];
    pos_ny_nz = [xv(i), yv(i) - eps, zv(i) - eps];
    
    BG_py_pz = cmag.FieldAndGradient(pos_py_pz, currents);
    BG_py_nz = cmag.FieldAndGradient(pos_py_nz, currents);
    BG_ny_pz = cmag.FieldAndGradient(pos_ny_pz, currents);
    BG_ny_nz = cmag.FieldAndGradient(pos_ny_nz, currents);
    
    % dBx/dy
    Bx_y(:) = (BG_py(1) - BG_ny(1))/(2);
    
    % dBx/dz
    Bx_z(:) = (BG_pz(1) - BG_nz(1))/2;
    
    % dBy/dz
    By_z(:) = (BG_pz(2) - BG_nz(2))/2;
    
    % d2Bx/dydz
    Bx_yz(:,i) = (BG_py_pz(1) - BG_ny_pz(1) - BG_py_nz(1) + BG_ny_nz(1)) / (4);
end

D = zeros(64,1);
D(9:16)     = fields(:,1);
D(17:24)    = fields(:,2);
D(25:32)    = fields(:,3);
D(33:40)    = Bx_y(:,2);
D(41:48)    = Bx_z(:,3);
D(49:56)    = By_z(:,5);
D(57:64)    = Bx_yz;

%D = reshape(fields', 1, [])';

a_sol = M \ D;
a_sol = [0; a_sol];
A_sol = reshape(a_sol, [4 4 4]);

B_fun = symfun(-gradient(tricubic(A_sol,x,y,z), [x y z]), [x y z]);

% test in plane
N = 8;

xvt = linspace(0,1,N);
yvt = linspace(0,1,N);

[xvtg, yvtg] = meshgrid(xvt, yvt);

fields_interp = zeros(N, N, 3);
fields_valid = zeros(N, N, 3);
scalar_pot = zeros(N,N);

for i = 1:N
    for j = 1:N
        pos = SCALE * [xvtg(i,j), yvtg(i,j), 0];
        scalar_pot(i,j) = double(tricubic(A_sol, xvtg(i,j), yvtg(i,j), 0));
        fields_interp(i, j, :) = B_fun(xvtg(i, j), yvtg(i,j), 0);
        BG = cmag.FieldAndGradient(pos, currents);
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
%             BG = cmag.FieldAndGradient(pos, currents);
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
        
