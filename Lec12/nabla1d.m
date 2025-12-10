function df=nabla1d(f,a, b,Nd);
%Determination of derivative of periodic function f(x). Domain is [a b].
%f: discrete codomain values of the function f.
%Sampling is very important: x=( (0:2*N-1)/(2*N) )*(b-a) +a;
%[a b]: The beginning and the end of the domain interval of f.
%Nd: Nd order of the derivative

N=length(f)/2;  % Half of Number of sample points.
rate=abs(b-a);  % A rate number if the interval length is not equal to 1.

f_fourier=fftshift(fft(f));  % Fourier tranformation f in [0 1].
 % Shift zero-frequency component to center of spectrum.

 %Derivation in frequency space.
frequency=-N:N-1;
derivation(frequency+N+1)=(2*pi*i*frequency/rate); % (2*pi*i*k_j)
f_fourier=f_fourier.*derivation.^Nd;

%Derivative in coordinate space.
df=ifft(ifftshift(f_fourier));