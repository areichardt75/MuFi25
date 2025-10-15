% Mufi 2025 - Lecture 6. - Nonlinear equations
% Solving nonlinear equations / chemical reaction 
% 
% x^2+2*x*y = 4; y^2-x*y = 1

fun = @(x) [x(1).^2+2*x(1).*x(2)-4; x(2).^2-x(1).*x(2)-1];

x0 = [0;1];
opts=  optimoptions('fsolve','Display','iter');
[x,fval] = fsolve(fun, x0, opts);
