% An electrostatic problem 
% Constrained Electrostatic Nonlinear Optimization Using Optimization
% Variables
% 

% Problem Geometry
[X,Y] = meshgrid(-1:.01:1);
z1 = -abs(X) - abs(Y);
z2 = -1 - sqrt(1-X.^2-Y.^2);
z2 = real(z2);
W1 = z1; W2 = z2;
W2(z1<z2) = nan;
W1(z1<z2) = nan;
hand = figure;
set(gcf,'Color','w');
surf(X,Y,W2,'LineStyle','none');
hold on;
surf(X,Y,W1,'LineStyle','none');
view(-44, 18)

%% Define Problem Variables
