% Circuit theory problem in case sinusoidal excitation
% Us = Ucs*cos(om*t+fi) 
% at om = om0 :
% Z1 = -10*j; Z2 = 20*j; Z3 = -10*j;
% R1 = 5; R2 = 5;
clearvars;
Z1 = -10*j; Z2 = 20*j; Z3 = -10*j;R1 = 5; R2 = 5;
Us = 10;
om = 1;
%% A. 
% Solve at om=om0 frequency
A = [1/R1+1/Z1+1/Z2 -1/Z2;-1/Z2 1/Z2+1/Z3+1/R2];
B = [Us/R1; 0];
x = A \ B;
Ur = x(2); % response 

%% Plot response and excitation in real time
% x(t) = Re{ Xcs*exp(j*om*t)}
T = 2*pi/om;
time = 0:T/100:2*T;
ust = real(Us*cos(om*time));
uresp = real(Ur*cos(om*time));

figure; 
  plot(time, ust, 'k-','LineWidth',2);
  hold on;
  plot(time, uresp, 'r-','LineWidth',2);
  xlabel('time'); ylabel('voltage');
  


