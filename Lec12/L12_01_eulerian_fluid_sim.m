% 2D Eulerian Fluid Simulation
% Solves incompressible Navier-Stokes equations on a fixed grid

clearvars; close all; %clc;
tic;
% Simulation parameters
N = 2*64;              % Grid resolution
dt = 0.01;           % Time step
visc = 0.0001;       % Viscosity
diff = 0.0001;       % Diffusion rate
iterations = 20;     % Pressure solver iterations
total_steps = 500;   % Total simulation steps

% Initialize grids
u = zeros(N, N);     % x-velocity
v = zeros(N, N);     % y-velocity
u_prev = zeros(N, N);
v_prev = zeros(N, N);

dens = zeros(N, N);  % Density/dye field
dens_prev = zeros(N, N);

% Variables to save for postprocessing
dens_save = zeros(total_steps, N, N);
u_save = zeros(total_steps, N, N);
v_save = zeros(total_steps, N, N);

% Add initial velocity disturbance
u(N/2-5:N/2+5, N/2-5:N/2+5) = 5;
v(N/2-5:N/2+5, N/2-5:N/2+5) = 3;

% Add initial density/dye
dens(N/2-3:N/2+3, N/2-3:N/2+3) = 1;

% Visualization setup
figure('Position', [100, 100, 800, 600]);

time1 = toc;
% Main simulation loop
for step = 1:total_steps
    % Add forces (example: rotating force field)
    [X, Y] = meshgrid(1:N, 1:N);
    cx = N/2; cy = N/2;
    force_x = 0.5 * (Y - cy) / N;
    force_y = -0.5 * (X - cx) / N;
    u = u + dt * force_x;
    v = v + dt * force_y;
    
    % Add density source
    if mod(step, 20) == 0
        dens(N/2-2:N/2+2, N/2-2:N/2+2) = dens(N/2-2:N/2+2, N/2-2:N/2+2) + 0.5;
    end
    
    % Velocity step
    [u, v] = velocity_step(u, v, u_prev, v_prev, visc, dt, N, iterations);
    
    % Density step
    dens = density_step(dens, dens_prev, u, v, diff, dt, N);
    
    % Visualization every 5 steps
    if mod(step, 2) == 0
        clf;

        % Plot density field
        subplot(1,2,1);
        imagesc(dens');
        colormap('hot');
        colorbar;
        title(sprintf('Density Field (Step %d)', step));
        axis equal tight;

        % Plot velocity field
        subplot(1,2,2);
        skip = 4;
        quiver(X(1:skip:end, 1:skip:end), Y(1:skip:end, 1:skip:end), ...
               u(1:skip:end, 1:skip:end)', v(1:skip:end, 1:skip:end)', 2);
        title('Velocity Field');
        axis equal tight;
        xlim([1 N]); ylim([1 N]);

        drawnow;
    end
    
    % Store previous values
    u_prev = u;
    v_prev = v;
    dens_prev = dens;
    u_save(step, :,:) = u;
    v_save(step, :,:) = v;
    dens_save(step, :,:) = dens;
end
time2 = toc;
fprintf('Time to initialize : %8.3f\nTime to solve : %8.3f\n',time1, time2);

%% Additional Functions

function [u, v] = velocity_step(u, v, u0, v0, visc, dt, N, iter)
    % Velocity step: diffusion, advection, and projection
    
    % Add forces
    u = u + dt * u0;
    v = v + dt * v0;
    
    % Diffusion
    u = diffuse(u, visc, dt, N, 1);
    v = diffuse(v, visc, dt, N, 2);
    
    % Project to make divergence-free
    [u, v] = project(u, v, N, iter);
    
    % Advection (self-advection)
    u = advect(u, u, v, dt, N, 1);
    v = advect(v, u, v, dt, N, 2);
    
    % Project again
    [u, v] = project(u, v, N, iter);
end

function dens = density_step(dens, dens0, u, v, diff, dt, N)
    % Density step: diffusion and advection
    
    % Add sources
    dens = dens + dt * dens0;
    
    % Diffusion
    dens = diffuse(dens, diff, dt, N, 0);
    
    % Advection
    dens = advect(dens, u, v, dt, N, 0);
end

function x = diffuse(x, diff, dt, N, b)
    % Implicit diffusion solver
    a = dt * diff * (N-2) * (N-2);
    x = lin_solve(x, x, a, 1 + 4*a, N, b, 20);
end

function x = advect(d, u, v, dt, N, b)
    % Semi-Lagrangian advection
    dt0 = dt * (N-2);
    x = d;
    
    for i = 2:N-1
        for j = 2:N-1
            % Backtrace
            x_pos = i - dt0 * u(i,j);
            y_pos = j - dt0 * v(i,j);
            
            % Clamp to grid
            x_pos = max(1.5, min(N-1.5, x_pos));
            y_pos = max(1.5, min(N-1.5, y_pos));
            
            % Bilinear interpolation
            i0 = floor(x_pos); i1 = i0 + 1;
            j0 = floor(y_pos); j1 = j0 + 1;
            
            s1 = x_pos - i0; s0 = 1 - s1;
            t1 = y_pos - j0; t0 = 1 - t1;
            
            x(i,j) = s0 * (t0 * d(i0,j0) + t1 * d(i0,j1)) + ...
                     s1 * (t0 * d(i1,j0) + t1 * d(i1,j1));
        end
    end
    
    x = set_bnd(x, N, b);
end

function [u, v] = project(u, v, N, iter)
    % Projection step to enforce incompressibility
    div = zeros(N, N);
    p = zeros(N, N);
    h = 1.0 / (N-2);
    
    % Compute divergence
    for i = 2:N-1
        for j = 2:N-1
            div(i,j) = -0.5 * h * (u(i+1,j) - u(i-1,j) + v(i,j+1) - v(i,j-1));
        end
    end
    div = set_bnd(div, N, 0);
    
    % Solve for pressure
    p = lin_solve(p, div, 1, 4, N, 0, iter);
    
    % Subtract pressure gradient
    for i = 2:N-1
        for j = 2:N-1
            u(i,j) = u(i,j) - 0.5 * (p(i+1,j) - p(i-1,j)) / h;
            v(i,j) = v(i,j) - 0.5 * (p(i,j+1) - p(i,j-1)) / h;
        end
    end
    u = set_bnd(u, N, 1);
    v = set_bnd(v, N, 2);
end

function x = lin_solve(x, x0, a, c, N, b, iter)
    % Gauss-Seidel linear solver
    for k = 1:iter
        for i = 2:N-1
            for j = 2:N-1
                x(i,j) = (x0(i,j) + a * (x(i-1,j) + x(i+1,j) + ...
                         x(i,j-1) + x(i,j+1))) / c;
            end
        end
        x = set_bnd(x, N, b);
    end
end

function x = set_bnd(x, N, b)
    % Set boundary conditions
    % b = 0: continuity, b = 1: horizontal velocity, b = 2: vertical velocity
    
    for i = 2:N-1
      if b==1 
        x(1,i) = -x(2,i);
        x(N,i) = -x(N-1,i);
      else
        x(1,i) = x(2,i); 
        x(N,i) = x(N-1,i);
      end
      if b==2
        x(i,1) = -x(i,2); 
        x(i,N) = -x(i,N-1);
      else 
        x(i,1) = x(i,2);
        x(i,N) = x(i,N-1);
      end

    end
    
    % Corners
    x(1, 1) = 0.5 * (x(2, 1) + x(1, 2));
    x(1, N) = 0.5 * (x(2, N) + x(1, N-1));
    x(N, 1) = 0.5 * (x(N-1, 1) + x(N, 2));
    x(N, N) = 0.5 * (x(N-1, N) + x(N, N-1));
end