% Mufi 2025 - Lecture 6. - Nonlinear equations
% Blackbody radiation
% 
clearvars;
clc;

fun = @(x) 5*(1-exp(-x))-x;
x0 = 2;
opts = optimoptions('fsolve','Display','iter-detailed');
xsol = fsolve(fun, x0,opts);
