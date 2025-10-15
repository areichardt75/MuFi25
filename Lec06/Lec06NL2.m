% Mufi 2025 - Lecture 6. - Nonlinear equations
% Solution of a nonlinear circuit (diode)
% un[V] = f(in), in[mA] 
% un = in^3-5*in^2+6*in
% x --> in 
% funchar(x) = x^3-5*x^2+6*x

% define function of characteristics
funchar = @(x) x.^3-5*x.^2+6*x;
% equation to solve for 
fun = @(x) (2/5)*x.^3-2*x.^2+17/5*x-2;
opts = optimoptions('fsolve','Display','iter');
x0 = 0;
xsol = fsolve(fun, x0, opts);

%% Visualization
in = 0:0.01:5;
un = funchar(in);
figure; 
  plot(in, un, 'r-');
  xlabel('in'); ylabel('un');
  title('Characteristics of diode');

hold on;
  plot(xsol, funchar(xsol), 'g.','MarkerSize',15);

hold off;



