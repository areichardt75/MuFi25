
clearvars; clc;

fprintf('Starting... ');

filename_to_load = 'Lagrangian1.mat';

load(filename_to_load);
[~,~,stepnumber] = size(densities_save);

fprintf('number of timesteps : %5d\n',stepnumber);

figure; 

