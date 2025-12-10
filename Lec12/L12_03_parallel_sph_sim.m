% Parallel SPH (Smoothed Particle Hydrodynamics) Fluid Simulation
% Uses MATLAB Parallel Computing Toolbox for multi-core acceleration

clearvars; close all; clc;
fname = 'SPH_parallel_dataX.mat';   % name of the file to store
%% Check for Parallel Computing Toolbox
if ~license('test', 'Distrib_Computing_Toolbox')
    warning('Parallel Computing Toolbox not available. Running in serial mode.');
    use_parallel = false;
else
    % Start parallel pool if not already running
    if isempty(gcp('nocreate'))
        parpool('local'); % Start parallel pool with default workers
    end
    use_parallel = true;
    fprintf('Using %d parallel workers\n', gcp().NumWorkers);
end

%% Simulation Parameters
num_particles = 2000;        % More particles to show parallel benefit
dt = 0.0008;                 % Time step
total_steps = 1000;          % Total simulation steps
h = 0.04;                    % Smoothing kernel radius
rest_density = 1000;         % Rest density (kg/m^3)
gas_constant = 2000;         % Pressure constant
viscosity = 250;             % Viscosity coefficient
particle_mass = 0.02;        % Mass per particle
gravity = [0, -9.8];         % Gravity acceleration

% Domain boundaries
x_min = 0; x_max = 1;
y_min = 0; y_max = 1;
damping = 0.5;

%% Initialize Particles
particles_per_row = round(sqrt(num_particles));
spacing = 0.015;
positions = zeros(num_particles, 2);
velocities = zeros(num_particles, 2);
positions_to_store = zeros(num_particles,2,total_steps);
velocities_to_store = zeros(num_particles,2,total_steps);

% Initialize particle positions (dam break scenario)
idx = 1;
for i = 1:particles_per_row
    for j = 1:particles_per_row
        if idx <= num_particles
            positions(idx, :) = [0.1 + i*spacing, 0.2 + j*spacing];
            idx = idx + 1;
        end
    end
end

%% Precompute kernel constants
h2 = h * h;
h9 = h^9;
poly6_const = 315 / (64 * pi * h9);
spiky_grad_const = -45 / (pi * h^6);
visc_lap_const = 45 / (pi * h^6);

%% Spatial Hashing Setup (for efficient neighbor search)
cell_size = h;
grid_width = ceil((x_max - x_min) / cell_size);
grid_height = ceil((y_max - y_min) / cell_size);

%% Visualization Setup
figure('Position', [100, 100, 900, 800]);

%% Performance Tracking
compute_times = zeros(total_steps, 1);

%% Main Simulation Loop
tic;
for step = 1:total_steps
    step_start = tic;
    
    %% Spatial Hashing for Fast Neighbor Search
    % Assign particles to grid cells
    grid_indices = floor((positions - [x_min, y_min]) / cell_size) + 1;
    grid_indices = max(1, min(grid_indices, [grid_width, grid_height]));
    
    % Create cell arrays for spatial hashing
    cell_hash = grid_indices(:,1) + (grid_indices(:,2)-1) * grid_width;
    [sorted_hash, sort_idx] = sort(cell_hash);
    sorted_positions = positions(sort_idx, :);
    sorted_velocities = velocities(sort_idx, :);
    
    % Find cell boundaries
    cell_start = zeros(grid_width * grid_height, 1);
    cell_end = zeros(grid_width * grid_height, 1);
    for i = 1:num_particles
        if i == 1 || sorted_hash(i) ~= sorted_hash(i-1)
            cell_start(sorted_hash(i)) = i;
        end
        cell_end(sorted_hash(i)) = i;
    end
    
    %% 1. PARALLEL Density Computation
    densities = zeros(num_particles, 1);
    
    if use_parallel
        % Parallel loop over particles
        parfor i = 1:num_particles
            pos_i = positions(i, :);
            dens = 0;
            
            % Get cell coordinates
            cell_x = grid_indices(i, 1);
            cell_y = grid_indices(i, 2);
            
            % Search in 3x3 neighborhood of cells
            for dx = -1:1
                for dy = -1:1
                    nx = cell_x + dx;
                    ny = cell_y + dy;
                    
                    if nx >= 1 && nx <= grid_width && ny >= 1 && ny <= grid_height
                        cell_idx = nx + (ny-1) * grid_width;
                        
                        % Check all particles in this cell
                        for j = cell_start(cell_idx):cell_end(cell_idx)
                            if j == 0
                                break;
                            end
                            
                            r_vec = pos_i - sorted_positions(j, :);
                            r2 = sum(r_vec.^2);
                            
                            if r2 < h2
                                dens = dens + particle_mass * poly6_const * (h2 - r2)^3;
                            end
                        end
                    end
                end
            end
            densities(i) = dens;
        end
    else
        % Serial version
        for i = 1:num_particles
            for j = 1:num_particles
                r_vec = positions(i, :) - positions(j, :);
                r2 = sum(r_vec.^2);
                if r2 < h2
                    densities(i) = densities(i) + particle_mass * poly6_const * (h2 - r2)^3;
                end
            end
        end
    end
    
    % Compute pressures
    pressures = gas_constant * (densities - rest_density);
    
    %% 2. PARALLEL Force Computation
    forces = zeros(num_particles, 2);
    
    if use_parallel
        parfor i = 1:num_particles
            pos_i = positions(i, :);
            vel_i = velocities(i, :);
            pressure_force = [0, 0];
            viscosity_force = [0, 0];
            
            % Get cell coordinates
            cell_x = grid_indices(i, 1);
            cell_y = grid_indices(i, 2);
            
            % Search in 3x3 neighborhood
            for dx = -1:1
                for dy = -1:1
                    nx = cell_x + dx;
                    ny = cell_y + dy;
                    
                    if nx >= 1 && nx <= grid_width && ny >= 1 && ny <= grid_height
                        cell_idx = nx + (ny-1) * grid_width;
                        
                        for j = cell_start(cell_idx):cell_end(cell_idx)
                            if j == 0
                                break;
                            end
                            
                            orig_j = sort_idx(j);
                            if i == orig_j
                                continue;
                            end
                            
                            r_vec = pos_i - sorted_positions(j, :);
                            r = norm(r_vec);
                            
                            if r < h && r > 1e-6
                                % Pressure force
                                pressure_force = pressure_force - particle_mass * ...
                                    (pressures(i) + pressures(orig_j)) / (2 * densities(orig_j)) * ...
                                    spiky_grad_const * (h - r)^2 * (r_vec / r);
                                
                                % Viscosity force
                                vel_diff = sorted_velocities(j, :) - vel_i;
                                viscosity_force = viscosity_force + particle_mass * ...
                                    vel_diff / densities(orig_j) * visc_lap_const * (h - r);
                            end
                        end
                    end
                end
            end
            
            forces(i, :) = pressure_force + viscosity * viscosity_force + densities(i) * gravity;
        end
    else
        % Serial version
        for i = 1:num_particles
            pressure_force = [0, 0];
            viscosity_force = [0, 0];
            
            for j = 1:num_particles
                if i == j, continue; end
                
                r_vec = positions(i, :) - positions(j, :);
                r = norm(r_vec);
                
                if r < h && r > 1e-6
                    pressure_force = pressure_force - particle_mass * ...
                        (pressures(i) + pressures(j)) / (2 * densities(j)) * ...
                        spiky_grad_const * (h - r)^2 * (r_vec / r);
                    
                    vel_diff = velocities(j, :) - velocities(i, :);
                    viscosity_force = viscosity_force + particle_mass * ...
                        vel_diff / densities(j) * visc_lap_const * (h - r);
                end
            end
            
            forces(i, :) = pressure_force + viscosity * viscosity_force + densities(i) * gravity;
        end
    end
    
    %% 3. PARALLEL Integration
    if use_parallel
        parfor i = 1:num_particles
            acc = forces(i, :) / densities(i);
            velocities(i, :) = velocities(i, :) + dt * acc;
            positions(i, :) = positions(i, :) + dt * velocities(i, :);
        end
    else
        accelerations = forces ./ densities;
        velocities = velocities + dt * accelerations;
        positions = positions + dt * velocities;
    end
    
    %% 4. Boundary Conditions (vectorized)
    % X boundaries
    left_wall = positions(:, 1) < x_min;
    positions(left_wall, 1) = x_min;
    velocities(left_wall, 1) = -damping * velocities(left_wall, 1);
    
    right_wall = positions(:, 1) > x_max;
    positions(right_wall, 1) = x_max;
    velocities(right_wall, 1) = -damping * velocities(right_wall, 1);
    
    % Y boundaries
    bottom_wall = positions(:, 2) < y_min;
    positions(bottom_wall, 2) = y_min;
    velocities(bottom_wall, 2) = -damping * velocities(bottom_wall, 2);
    
    top_wall = positions(:, 2) > y_max;
    positions(top_wall, 2) = y_max;
    velocities(top_wall, 2) = -damping * velocities(top_wall, 2);
    
    compute_times(step) = toc(step_start);
    
    % %% 5. Visualization (every 10 steps)
    % if mod(step, 100) == 0
    %     cla;
    % 
    %     vel_mag = sqrt(sum(velocities.^2, 2));
    %     scatter(positions(:, 1), positions(:, 2), 20, densities, 'filled');
    %     colormap('jet');
    %     cb = colorbar;
    %     ylabel(cb, 'Density (kg/m^3)');
    % 
    %     xlim([x_min, x_max]);
    %     ylim([y_min, y_max]);
    %     axis equal;
    %     rectangle('Position', [x_min, y_min, x_max-x_min, y_max-y_min], ...
    %               'EdgeColor', 'k', 'LineWidth', 2);
    % 
    %     title(sprintf('Parallel SPH - Step %d/%d (%.2f ms/step)', ...
    %           step, total_steps, mean(compute_times(max(1,step-50):step))*1000));
    %     xlabel('X Position (m)');
    %     ylabel('Y Position (m)');
    % 
    %     text(0.02, 0.95, sprintf('Particles: %d', num_particles), 'Units', 'normalized');
    %     text(0.02, 0.90, sprintf('Parallel: %s', string(use_parallel)), 'Units', 'normalized');
    % 
    %     drawnow;
    % end
    positions_to_store(:,:,step) = positions(:,:);
    velocities_to_store(:,:,step)= velocities(:,:);
end
save(fname, 'positions_to_store','velocities_to_store');
total_time = toc;

%% Performance Report
fprintf('\n=== Performance Report ===\n');
fprintf('Total simulation time: %.2f seconds\n', total_time);
fprintf('Average time per step: %.2f ms\n', mean(compute_times)*1000);
fprintf('Particles: %d\n', num_particles);
fprintf('Parallel mode: %s\n', string(use_parallel));
if use_parallel
    fprintf('Number of workers: %d\n', gcp().NumWorkers);
end
fprintf('Speedup potential: ~%dx with GPU implementation\n', 10);

disp('Simulation complete!');