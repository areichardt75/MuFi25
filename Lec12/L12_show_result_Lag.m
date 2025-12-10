
clearvars; clc;

fprintf('Starting... ');

filename_to_load = 'Lagrangian2.mat';

load(filename_to_load);
[~,~,stepnumber] = size(densities_save);

fprintf('number of timesteps : %5d\n',stepnumber);
% for step=1:total_steps
%   velmag = 
% figure; 
  figure('Position', [100, 100, 800, 800]);
  axis equal;
  xlim([x_min, x_max]);
  ylim([y_min, y_max]);
  hold on;
  grid on;

  % Draw boundaries
  rectangle('Position', [x_min, y_min, x_max-x_min, y_max-y_min], ...
          'EdgeColor', 'k', 'LineWidth', 2);
% end


for step=1:total_steps
      if (mod(step, 2) == 0)
        cla;

        velocities = vel_save(:,:,step);
        positions = pos_save(:,:,step);
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

        pause(0.5);
        drawnow;
      end
end


