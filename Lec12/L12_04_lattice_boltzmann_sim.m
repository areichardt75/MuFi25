% 2D Lattice Boltzmann Method (LBM) - D2Q9 Model
% Simulates fluid flow using kinetic theory on a fixed lattice

clear; close all; clc;

%% Simulation Parameters
Nx = 200;                    % Lattice width
Ny = 100;                    % Lattice height
total_steps = 5000;          % Simulation steps
omega = 1.7;                 % Relaxation parameter (1 < omega < 2)
                             % Related to viscosity: nu = (1/omega - 0.5)/3

% Flow parameters
u_max = 0.1;                 % Maximum velocity (must be << 1 for stability)
Re = 100;                    % Reynolds number

%% D2Q9 Lattice (9 velocities in 2D)
% Velocity directions:
%   6   2   5
%     \ | /
%   3 - 0 - 1
%     / | \
%   7   4   8

% Discrete velocities
ex = [0,  1,  0, -1,  0,  1, -1, -1,  1];
ey = [0,  0,  1,  0, -1,  1,  1, -1, -1];

% Lattice weights
w = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];

% Opposite directions (for bounce-back)
opp = [1, 4, 5, 2, 3, 8, 9, 6, 7];

%% Initialize Distribution Functions
% f(x, y, direction) = particle distribution
f = zeros(Ny, Nx, 9);
feq = zeros(Ny, Nx, 9);

% Initialize with equilibrium distribution (stationary fluid)
rho = ones(Ny, Nx);          % Density
ux = zeros(Ny, Nx);          % X-velocity
uy = zeros(Ny, Nx);          % Y-velocity

% Set initial velocity (Poiseuille flow profile)
for j = 1:Ny
    ux(j, :) = u_max * 4 * (j-1) * (Ny-j) / (Ny^2);
end

% Initialize distribution functions to equilibrium
for i = 1:9
    feq(:, :, i) = equilibrium(rho, ux, uy, i, ex, ey, w);
    f(:, :, i) = feq(:, :, i);
end

%% Create Obstacle (Cylinder)
obstacle = false(Ny, Nx);
cx = Nx/4; cy = Ny/2; radius = Ny/8;
[X, Y] = meshgrid(1:Nx, 1:Ny);
obstacle = ((X - cx).^2 + (Y - cy).^2) < radius^2;

%% Visualization Setup
figure('Position', [100, 100, 1200, 500]);

%% Main LBM Loop
for step = 1:total_steps
    
    %% 1. Compute Macroscopic Variables
    rho = sum(f, 3);
    ux = zeros(Ny, Nx);
    uy = zeros(Ny, Nx);
    
    for i = 1:9
        ux = ux + f(:, :, i) * ex(i);
        uy = uy + f(:, :, i) * ey(i);
    end
    
    ux = ux ./ rho;
    uy = uy ./ rho;
    
    %% 2. Apply Boundary Conditions
    
    % Inlet (left boundary) - constant velocity
    for j = 1:Ny
        u_inlet = u_max * 4 * (j-1) * (Ny-j) / (Ny^2);
        ux(j, 1) = u_inlet;
        uy(j, 1) = 0;
        rho(j, 1) = 1;
    end
    
    % Outlet (right boundary) - constant pressure/density
    rho(:, Nx) = 1;
    ux(:, Nx) = ux(:, Nx-1);
    uy(:, Nx) = uy(:, Nx-1);
    
    % Top and bottom walls - no-slip (zero velocity)
    ux(1, :) = 0;
    uy(1, :) = 0;
    ux(Ny, :) = 0;
    uy(Ny, :) = 0;
    
    % Obstacle - no-slip
    ux(obstacle) = 0;
    uy(obstacle) = 0;
    
    %% 3. Compute Equilibrium Distribution
    for i = 1:9
        feq(:, :, i) = equilibrium(rho, ux, uy, i, ex, ey, w);
    end
    
    %% 4. Collision Step (BGK approximation)
    f = f - omega * (f - feq);
    
    %% 5. Streaming Step (propagate distributions)
    f_new = f;
    for i = 1:9
        % Stream in direction i
        f_new(:, :, i) = circshift(f(:, :, i), [ey(i), ex(i)]);
    end
    f = f_new;
    
    %% 6. Bounce-back on Obstacles
    for i = 1:9
        % Reverse distributions on obstacle nodes
        temp = f(:, :, i);
        temp(obstacle) = f(obstacle, opp(i));
        f(:, :, i) = temp;
    end
    
    %% 7. Bounce-back on Walls
    % Bottom wall
    f(1, :, 2) = f(1, :, 4);
    f(1, :, 5) = f(1, :, 7);
    f(1, :, 6) = f(1, :, 8);
    
    % Top wall
    f(Ny, :, 4) = f(Ny, :, 2);
    f(Ny, :, 7) = f(Ny, :, 5);
    f(Ny, :, 8) = f(Ny, :, 6);
    
    %% 8. Visualization (every 50 steps)
    if mod(step, 50) == 0
        % Compute velocity magnitude and vorticity
        u_mag = sqrt(ux.^2 + uy.^2);
        
        % Compute vorticity (curl of velocity field)
        [dux_dy, ~] = gradient(ux);
        [~, duy_dx] = gradient(uy);
        vorticity = duy_dx - dux_dy;
        
        % Plot velocity magnitude
        subplot(1, 2, 1);
        contourf(u_mag, 20, 'LineStyle', 'none');
        colormap('jet');
        colorbar;
        hold on;
        contour(double(obstacle), [0.5, 0.5], 'k', 'LineWidth', 2);
        hold off;
        title(sprintf('Velocity Magnitude (Step %d)', step));
        xlabel('X'); ylabel('Y');
        axis equal tight;
        
        % Plot vorticity
        subplot(1, 2, 2);
        contourf(vorticity, 20, 'LineStyle', 'none');
        colormap('jet');
        colorbar;
        hold on;
        contour(double(obstacle), [0.5, 0.5], 'k', 'LineWidth', 2);
        hold off;
        title('Vorticity (Flow Rotation)');
        xlabel('X'); ylabel('Y');
        axis equal tight;
        
        drawnow;
    end
end

disp('Simulation complete!');

%% Helper Function: Equilibrium Distribution
function feq_i = equilibrium(rho, ux, uy, i, ex, ey, w)
    % Compute equilibrium distribution for direction i
    % This is the Maxwell-Boltzmann-like distribution
    
    eu = ex(i) * ux + ey(i) * uy;
    u2 = ux.^2 + uy.^2;
    
    feq_i = rho * w(i) .* (1 + 3*eu + 4.5*eu.^2 - 1.5*u2);
end