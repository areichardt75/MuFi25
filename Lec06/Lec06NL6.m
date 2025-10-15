% Mufi 2025 - Lecture 6. - Nonlinear equations
% Nonlinear circuit with two diodes
% u2 + i1*R1 + i2*(R1+R2) - U0 = 0
% u1 + i1*(R1+R3) + i2*R1 - U0 = 0
% u1(i1) = 2*i1^2-1
% u2(i2) = (2-i2)*(1+i2) = -i2^2+i2+2
% x --> u1, i1, u2, i2
diode1 = @(x) 2*x.^2-1;
diode2 = @(x) -x.^2+x+2;
U0 = 10;
R1 = 5; 
R2 = 10;
R3 = 10;
nlin = @(x) [x(3)+x(2)*R2+x(4)*(R1+R2)-U0;...
  x(1)+x(2)*(R1+R3)+x(4)*R1-U0;...
  x(1)-diode1(x(2));...
  x(3)-diode2(x(4))];

x0 = [0;0;0;0];
opts = optimoptions('fsolve','Display','iter');
x = fsolve(nlin, x0, opts)







function res = nlinfun(x)
U0 = 10;
R1 = 5; 
R2 = 10;
R3 = 10;
res = [x(3)+x(2)*R2+x(4)*(R1+R2)-U0;...
  x(1)+x(2)*(R1+R3)+x(4)*R1-U0;...
  x(1)-2*x(2).^2+1;...
  x(3)+x(4).^2-x(4)-2];
end


