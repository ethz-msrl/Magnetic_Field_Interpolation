% Inspired by "Tricubic Interpolation in three dimensions" from F. Lekien
% and J. Marsden

% The idea is to interpolate on a 3D cube by fixing continuity in 
% [f, df_dx, df_dy, df_dz, d2f_dxdy, d2f_dxdz, d2f_dydz, d3f_dxdydz]
% We interpolate one scalar field per component Bx, By, Bz

%% Interpolation part
% This param scales the size of the cube
SCALE = 30e-3;

% the cube corners
xv = SCALE * [0, 0, 0, 0, 1, 1, 1, 1];
yv = SCALE * [0, 0, 1, 1, 0, 0, 1, 1];
zv = SCALE * [0, 1, 0, 1, 0, 1, 0, 1];

% currents running from the MMS
%currents = 16 * rand(8,1) - 8;
currents = [6.7573    4.3353   -7.3174   -1.9490    3.2694    3.6722   -4.4116   -3.6951]';

fields = zeros(3, 8);

for i = 1:8
    pos = [xv(i), yv(i), zv(i)];
    BG = magSystem.FieldAndGradient(pos, currents);
    fields(:, i) = SCALE*BG(1:3);
end

% Here we calculate the function values and derivatives using finite
% differences

eps = SCALE;
B_x = zeros(3,8);
B_y = zeros(3,8);
B_z = zeros(3,8);
B_xy = zeros(3,8);
B_xz = zeros(3,8);
B_yz = zeros(3,8);
B_xyz = zeros(3,8);

for i = 1:8
    pos = [xv(i), yv(i), zv(i)];
    BG = magSystem.FieldAndGradient(pos, currents);
    pos_px = [xv(i) + eps, yv(i), zv(i)];
    pos_nx = [xv(i) - eps, yv(i), zv(i)];
    BG_px = magSystem.FieldAndGradient(pos_px, currents);
    BG_nx = magSystem.FieldAndGradient(pos_nx, currents);
    B_x(:,i) = (BG_px(1:3) - BG_nx(1:3))/(2);
    
    pos_py = [xv(i), yv(i) + eps, zv(i)];
    pos_ny = [xv(i), yv(i) - eps, zv(i)];
    BG_py = magSystem.FieldAndGradient(pos_py, currents);
    BG_ny = magSystem.FieldAndGradient(pos_ny, currents);
    B_y(:,i) = (BG_py(1:3) - BG_ny(1:3))/(2);
    
    pos_pz = [xv(i), yv(i), zv(i) + eps];
    pos_nz = [xv(i), yv(i), zv(i) - eps];
    BG_pz = magSystem.FieldAndGradient(pos_pz, currents);
    BG_nz = magSystem.FieldAndGradient(pos_nz, currents);
    B_z(:,i) = (BG_pz(1:3) - BG_nz(1:3))/(2);
    
    pos_px_py = [xv(i) + eps, yv(i) + eps, zv(i)];
    pos_px_ny = [xv(i) + eps, yv(i) - eps, zv(i)];
    pos_nx_py = [xv(i) - eps, yv(i) + eps, zv(i)];
    pos_nx_ny = [xv(i) - eps, yv(i) - eps, zv(i)];
    
    BG_px_py = magSystem.FieldAndGradient(pos_px_py, currents);
    BG_px_ny = magSystem.FieldAndGradient(pos_px_ny, currents);
    BG_nx_py = magSystem.FieldAndGradient(pos_nx_py, currents);
    BG_nx_ny = magSystem.FieldAndGradient(pos_nx_ny, currents);
    
    pos_px_pz = [xv(i) + eps, yv(i), zv(i) + eps];
    pos_px_nz = [xv(i) + eps, yv(i), zv(i) - eps];
    pos_nx_pz = [xv(i) - eps, yv(i), zv(i) + eps];
    pos_nx_nz = [xv(i) - eps, yv(i), zv(i) - eps];
    
    BG_px_pz = magSystem.FieldAndGradient(pos_px_pz, currents);
    BG_px_nz = magSystem.FieldAndGradient(pos_px_nz, currents);
    BG_nx_pz = magSystem.FieldAndGradient(pos_nx_pz, currents);
    BG_nx_nz = magSystem.FieldAndGradient(pos_nx_nz, currents);
    
    pos_py_pz = [xv(i), yv(i) + eps, zv(i) + eps];
    pos_py_nz = [xv(i), yv(i) + eps, zv(i) - eps];
    pos_ny_pz = [xv(i), yv(i) - eps, zv(i) + eps];
    pos_ny_nz = [xv(i), yv(i) - eps, zv(i) - eps];
    
    BG_py_pz = magSystem.FieldAndGradient(pos_py_pz, currents);
    BG_py_nz = magSystem.FieldAndGradient(pos_py_nz, currents);
    BG_ny_pz = magSystem.FieldAndGradient(pos_ny_pz, currents);
    BG_ny_nz = magSystem.FieldAndGradient(pos_ny_nz, currents);
    
    B_xy(:,i) = (BG_px_py(1:3) - BG_nx_py(1:3) - BG_px_ny(1:3) + BG_nx_ny(1:3)) / (4);
    B_xz(:,i) = (BG_px_pz(1:3) - BG_nx_pz(1:3) - BG_px_nz(1:3) + BG_nx_nz(1:3)) / (4);
    B_yz(:,i) = (BG_py_pz(1:3) - BG_ny_pz(1:3) - BG_py_nz(1:3) + BG_ny_nz(1:3)) / (4);
    
    pos_px_py_pz = [xv(i) + eps, yv(i) + eps, zv(i) + eps];
    pos_nx_py_pz = [xv(i) - eps, yv(i) + eps, zv(i) + eps];
    pos_px_ny_pz = [xv(i) + eps, yv(i) - eps, zv(i) + eps];
    pos_nx_ny_pz = [xv(i) - eps, yv(i) - eps, zv(i) + eps];
    pos_px_py_nz = [xv(i) + eps, yv(i) + eps, zv(i) - eps];
    pos_nx_py_nz = [xv(i) - eps, yv(i) + eps, zv(i) - eps];
    pos_px_ny_nz = [xv(i) + eps, yv(i) - eps, zv(i) - eps];
    pos_nx_ny_nz = [xv(i) - eps, yv(i) - eps, zv(i) - eps];
    
    BG_px_py_pz = magSystem.FieldAndGradient(pos_px_py_pz, currents);
    BG_nx_py_pz = magSystem.FieldAndGradient(pos_nx_py_pz, currents);
    BG_px_ny_pz = magSystem.FieldAndGradient(pos_px_ny_pz, currents);
    BG_nx_ny_pz = magSystem.FieldAndGradient(pos_nx_ny_pz, currents);
    BG_px_py_nz = magSystem.FieldAndGradient(pos_px_py_nz, currents);
    BG_nx_py_nz = magSystem.FieldAndGradient(pos_nx_py_nz, currents);
    BG_px_ny_nz = magSystem.FieldAndGradient(pos_px_ny_nz, currents);
    BG_nx_ny_nz = magSystem.FieldAndGradient(pos_nx_ny_nz, currents);
    
    B_xyz(:,i) = (BG_px_py_pz(1:3) - BG_nx_py_pz(1:3) - BG_px_ny_pz(1:3) ...
        + BG_nx_ny_pz(1:3) - BG_px_py_nz(1:3) - BG_nx_py_nz(1:3) + BG_px_ny_nz(1:3) ...
        - BG_nx_ny_nz(1:3))/(8); 
end

% interpolating x field
D = zeros(64,1);
D(1:8)      = fields(1,:);
D(9:16)     = B_x(1,:);
D(17:24)    = B_y(1,:);
D(25:32)    = B_z(1,:);
D(33:40)    = B_xy(1,:);
D(41:48)    = B_xz(1,:);
D(49:56)    = B_yz(1,:);
D(57:end)   = B_xyz(1,:);

a_sol = M \ D;
A_sol = reshape(a_sol, [4 4 4]);

Bx_fun = symfun(tricubic(A_sol,x,y,z), [x y z]);

% interpolating y field
D = zeros(64,1);
D(1:8)      = fields(2,:);
D(9:16)     = B_x(2,:);
D(17:24)    = B_y(2,:);
D(25:32)    = B_z(2,:);
D(33:40)    = B_xy(2,:);
D(41:48)    = B_xz(2,:);
D(49:56)    = B_yz(2,:);
D(57:end)   = B_xyz(2,:);

a_sol = M \ D;
A_sol = reshape(a_sol, [4 4 4]);

By_fun = symfun(tricubic(A_sol,x,y,z), [x y z]);

% interpolating z field
D = zeros(64,1);
D(1:8)      = fields(3,:);
D(9:16)     = B_x(3,:);
D(17:24)    = B_y(3,:);
D(25:32)    = B_z(3,:);
D(33:40)    = B_xy(3,:);
D(41:48)    = B_xz(3,:);
D(49:56)    = B_yz(3,:);
D(57:end)   = B_xyz(3,:);

a_sol = M \ D;
A_sol = reshape(a_sol, [4 4 4]);

Bz_fun = symfun(tricubic(A_sol,x,y,z), [x y z]);

% test in cube
N = 8;
xvt = linspace(0,1,N);
yvt = linspace(0,1,N);
zvt = linspace(0,1,N);

[xvtg, yvtg, zvtg] = meshgrid(xvt, yvt, zvt);

fields_interp = zeros(N, N, N, 3);
fields_valid = zeros(N, N, N, 3);

for i = 1:N
    for j = 1:N
        for k =1:N
            pos = SCALE* [xvt(i), yvt(j), zvt(k)];
            fields_interp(i, j, k, 1) = Bx_fun(xvt(i), yvt(j), zvt(k));
            fields_interp(i, j, k, 2) = By_fun(xvt(i), yvt(j), zvt(k));
            fields_interp(i, j, k, 3) = Bz_fun(xvt(i), yvt(j), zvt(k));
            BG = magSystem.FieldAndGradient(pos, currents);
            fields_valid(i,j,k,:) = BG(1:3);
        end
    end
end

plot_quiver_2d;
