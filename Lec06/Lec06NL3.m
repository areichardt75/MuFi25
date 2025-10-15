% Mufi 2025 - Lecture 6. - Nonlinear equations
% Solution of Van der Waals gas in case of oxygene gas
% 
P = 10; 
T = 300;
R = 0.08206; % in atm*L/(mol*K) 1atm = 101300 Pa!
a = 1.37;
b = 0.0318;

% definition of state equation in full form
fun1 = @(V) (P+a./(V.^2)).*(V-b)-R*T;

% definition of state equation in rearranged form
% x--> V
fun2 = @(x) P*x.^3-(P*b+R*T)*x.^2+a*x-a*b;

% starting point
x0 = -5;
% options
opts = optimoptions('fsolve', 'Display','iter');

% Start stopwatch for measuring time
% t0 = tic;
% let's solve it using original form
v1 = fsolve(fun1, x0, opts);
fprintf('Solution of fun1 : V = %.5f\n',v1);
% t2 = toc; %stop time at the end of 1st part

% and solve it in rearranged form
v2 = fsolve(fun2, x0, opts);
fprintf('Solution of fun2 : V = %.5f\n',v2);
% t3 = toc; % stop time at the end of 2nd part

% fprintf('Times :\n1. run : %.6f sec\n2.run : %.6f sec\n',t2-t0,t3-t2);