% Circuit theory problem in case of sinusoidal excitation
% Us = Ucs*cos(om*t+fi) 
% at om = om0 :
% Z1 = -10*j; Z2 = 20*j; Z3 = -10*j;
% R1 = 5; R2 = 5;
% 
% Problem 1/B. Frequency dependence of response
% om = om/100:...:om*100;
% xi = om/om0
xi = logspace(-2,2,1e4);
Uresp = zeros(size(xi));
for idx=1:length(xi)
  Ut = filterLC(xi(idx), 10);
  Uresp(idx) = Ut/10;
end

figure; plot(xi, abs(Uresp), 'r-');
xlabel('frequency multiplication');
ylabel('transfer amplitude');
title('Voltage amplification');

% Poor quality of graph! 

%% Use logarithmic axis on frequency axis
% 
figure;
  semilogx(xi, abs(Uresp), 'k-','LineWidth',2);
  xlabel('frequency multiplication [logarithmic]');
  ylabel('transfer amplitude');
  title('Voltage amplification');
  ax = gca;
  ax.XGrid = 'on'; ax.YGrid = 'on';
  



