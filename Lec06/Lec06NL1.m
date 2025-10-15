% Mufi 2025 - Lecture 6. - Nonlinear equations
% Solving set of nonlinear equations 
% stored in pfun(x,c) equation
% We solve equations for c=-1
c = -1;
fun = @(x) pfun(x,c);
% solve system starting at [0,1]
x0 = [0;1];
x0 = [-2;2];
% set options to see what is happening inside fsolve
% 
opts = optimoptions('fsolve','Display','iter');
x = fsolve(fun, x0, opts);










function F = pfun(x,c)
  F = [2*x(1)+x(2)-exp(c*x(1));-x(1)+2*x(2)-exp(c*x(2))];
end
