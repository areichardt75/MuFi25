% 2D SPH (Smoothed Particle Hydrodynamics) Lagrangian Fluid Simulation
% Tracks individual fluid particles and computes interactions

clear; close all; clc;

visualize = 10;
%% Simulation Parameters
num_particles = 500;        % Number of fluid particles
dt = 0.001;                  % Time step
total_steps = 2000;          % Total simulation steps
h = 0.05;                    % Smoothing kernel radius
rest_density = 1000;         % Rest density (kg/m^3)
gas_constant = 2000;         % Pressure constant
viscosity = 250;             % Viscosity coefficient
particle_mass = 0.02;        % Mass per particle
gravity = [0, -9.8];         % Gravity acceleration

% Domain boundaries
x_min = 0; x_max = 1;
y_min = 0; y_max = 1;
damping = 0.5;               % Boundary damping coefficient

%% Initialize Particles
% Create particles in a block
particles_per_row = round(sqrt(num_particles));
spacing = 0.02;
positions = zeros(num_particles, 2);
velocities = zeros(num_particles, 2);
forces = zeros(num_particles, 2);
densities = zeros(num_particles, 1);
pressures = zeros(num_particles, 1);
pos_save = zeros(num_particles,2, total_steps);
vel_save = zeros(num_particles,2, total_steps);
forces_save = zeros(num_particles,2, total_steps);
densities_save = zeros(num_particles,1, total_steps);
pressures_save = zeros(num_particles,1, total_steps);


% Initialize particle positions (water block)
idx = 1;
for i = 1:particles_per_row
    for j = 1:particles_per_row
        if idx <= num_particles
            positions(idx, :) = [0.2 + i*spacing, 0.5 + j*spacing];
            idx = idx + 1;
        end
    end
end

% Add some initial velocity (optional)
velocities(:, 1) = 2;  % Initial horizontal velocity

%% Precompute kernel constants
h2 = h * h;
h9 = h^9;
poly6_const = 315 / (64 * pi * h9);
spiky_grad_const = -45 / (pi * h^6);
visc_lap_const = 45 / (pi * h^6);

%% Visualization Setup
if visualize~=0
  figure('Position', [100, 100, 800, 800]);
  axis equal;
  xlim([x_min, x_max]);
  ylim([y_min, y_max]);
  hold on;
  grid on;

  % Draw boundaries
  rectangle('Position', [x_min, y_min, x_max-x_min, y_max-y_min], ...
          'EdgeColor', 'k', 'LineWidth', 2);
end
%% Main Simulation Loop
for step = 1:total_steps
    
    %% 1. Compute Densities and Pressures
    densities = zeros(num_particles, 1);
    
    for i = 1:num_particles
        for j = 1:num_particles
            r_vec = positions(i, :) - positions(j, :);
            r2 = sum(r_vec.^2);
            
            if r2 < h2
                % Poly6 kernel for density
                densities(i) = densities(i) + particle_mass * ...
                    poly6_const * (h2 - r2)^3;
            end
        end
    end
    
    % Compute pressures using equation of state
    pressures = gas_constant * (densities - rest_density);
    
    %% 2. Compute Forces
    forces = zeros(num_particles, 2);
    
    for i = 1:num_particles
        pressure_force = [0, 0];
        viscosity_force = [0, 0];
        
        for j = 1:num_particles
            if i == j
                continue;
            end
            
            r_vec = positions(i, :) - positions(j, :);
            r = norm(r_vec);
            
            if r < h && r > 1e-6
                % Pressure force (Spiky kernel gradient)
                pressure_force = pressure_force - particle_mass * ...
                    (pressures(i) + pressures(j)) / (2 * densities(j)) * ...
                    spiky_grad_const * (h - r)^2 * (r_vec / r);
                
                % Viscosity force (Laplacian)
                vel_diff = velocities(j, :) - velocities(i, :);
                viscosity_force = viscosity_force + particle_mass * ...
                    vel_diff / densities(j) * visc_lap_const * (h - r);
            end
        end
        
        % Total force = pressure + viscosity + gravity
        forces(i, :) = pressure_force + viscosity * viscosity_force + ...
                       densities(i) * gravity;
    end
    
    %% 3. Integrate (Explicit Euler)
    accelerations = forces ./ densities;
    velocities = velocities + dt * accelerations;
    positions = positions + dt * velocities;
    
    %% 4. Boundary Conditions
    for i = 1:num_particles
        % X boundaries
        if positions(i, 1) < x_min
            positions(i, 1) = x_min;
            velocities(i, 1) = -damping * velocities(i, 1);
        elseif positions(i, 1) > x_max
            positions(i, 1) = x_max;
            velocities(i, 1) = -damping * velocities(i, 1);
        end
        
        % Y boundaries
        if positions(i, 2) < y_min
            positions(i, 2) = y_min;
            velocities(i, 2) = -damping * velocities(i, 2);
        elseif positions(i, 2) > y_max
            positions(i, 2) = y_max;
            velocities(i, 2) = -damping * velocities(i, 2);
        end
    end
    
    %% 5. Visualization (every 10 steps)
    % Now, we are focused on computing solution of problem and save results
    % to postprocess later. 
    % 
    if (visualize ~=0) && (mod(step, 2) == 0)
        cla;

        % Color particles by velocity magnitude
        vel_mag = sqrt(sum(velocities.^2, 2));
        scatter(positions(:, 1), positions(:, 2), 30, vel_mag, 'filled');
        colormap('jet');
        colorbar;
        caxis([0, max(vel_mag)+0.1]);

        % Redraw boundaries
        rectangle('Position', [x_min, y_min, x_max-x_min, y_max-y_min], ...
                  'EdgeColor', 'k', 'LineWidth', 2);

        title(sprintf('SPH Simulation - Step %d/%d', step, total_steps));
        xlabel('X Position (m)');
        ylabel('Y Position (m)');

        % Display stats
        text(0.02, 0.95, sprintf('Particles: %d', num_particles), ...
             'Units', 'normalized');
        text(0.02, 0.90, sprintf('Avg Density: %.1f kg/m^3', mean(densities)), ...
             'Units', 'normalized');
        text(0.02, 0.85, sprintf('Max Velocity: %.2f m/s', max(vel_mag)), ...
             'Units', 'normalized');

        drawnow;
    end
    %% 6. Save results to variables
      pos_save(:,:,step) = positions(:,:);
      vel_save(:,:,step) = velocities(:,:); 
      forces_save(:,:,step) = forces(:,:); 
      densities_save(:,:,step) = densities(:,:);
      pressures_save(:,:,step) = pressures(:,:);

    if mod(step,10) == 0 
      fprintf('%5d | ', step);
      if mod(step,100) == 0 
        fprintf('\n');
      end 
    end
end

disp('Simulation complete!');

%% Additional Helper Functions for Analysis

function W = poly6_kernel(r, h)
    % Poly6 smoothing kernel
    if r >= 0 && r <= h
        const = 315 / (64 * pi * h^9);
        W = const * (h^2 - r^2)^3;
    else
        W = 0;
    end
end

function grad_W = spiky_gradient(r_vec, h)
    % Spiky kernel gradient (for pressure)
    r = norm(r_vec);
    if r >= 0 && r <= h && r > 1e-6
        const = -45 / (pi * h^6);
        grad_W = const * (h - r)^2 * (r_vec / r);
    else
        grad_W = [0, 0];
    end
end

function lap_W = viscosity_laplacian(r, h)
    % Viscosity kernel Laplacian
    if r >= 0 && r <= h
        const = 45 / (pi * h^6);
        lap_W = const * (h - r);
    else
        lap_W = 0;
    end
end