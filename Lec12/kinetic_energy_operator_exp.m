function KO=kinetic_energy_operator_exp(f,a, b, m, dt)
%Determination of exponetial function of kinetic energy operator.
%f: discrete codomain values of the function f.
%Sampling is very important: x=( (0:2*N-1)/(2*N) )*(b-a) + a;
%[a b]: The beginning and the end of the domain interval of f.
%Nd: if (Nd >=0) order of the derivative
% if (Nd<0) integral of f

N=length(f)/2; % Half of Number of sample points.
rate=abs(b-a); % A rate number if the interval length is not equal to 1.

f_fourier=fftshift(fft(f)); % Fourier tranformation f in [0 1].
 % Shift zero-frequency component to center of spectrum.

%Expontential function of Kinetic energy operator in momentum space.
frequency=-N:N-1;
derivation(frequency+N+1)=(2*pi*i*frequency/rate); % (2*pi*i*k_j)
f_fourier= exp ( 0.5*sqrt(-1)*dt/2/m * derivation.^2 ).* f_fourier;


%Derivative in coordinate space.
KO=ifft(ifftshift(f_fourier));