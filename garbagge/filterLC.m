function [Uresp] = filterLC(xi, Us)
%FILTERLC 
%   
Z10 = -10*j; 
Z20 = 20*j; 
Z30 = -10*j; 
R1 = 5;
R2 = 5;
Z1 = Z10/xi; Z2 = Z20*xi; Z3 = Z30/xi;
A = [1/R1+1/Z1+1/Z2 -1/Z2;-1/Z2 1/Z2+1/Z3+1/R2];
B = [Us/R1; 0];
x = A\B;
Uresp = x(2);
end

