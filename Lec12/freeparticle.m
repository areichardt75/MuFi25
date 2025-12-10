%Free particle quantum mechanical motion.
%If wave number k=0, analytical and numerical solutions are compared.
%Free particle ensemble is described by Gaussina wave-packet.
%
%Solution of timedependent Scrodinger equation by Split Operator Method (SPO).
%exp(-i/h_bar H dt) ~ exp(-i/h_bar/2 K dt) exp (-h_bar V dt) exp(-i/h_bar/2 K dt)
%ih dPsi/dt = H Psi
%H = -h_bar^2/2/m d^2/dx^2 + V
%
%References:
%C. Leforester and coworkers, Journal of computational physics, vol 94, 59-80 (1991)
%G. Varga, Journal of Physics: Condensed Matter, vol 14 (2002) 6081-6107
%
%%Atomic units are used:
%	Planck constant: h=6.6252e-34 Js, h_bar=1.
%	Mass: mu=9.109558e-31 kg; Electron mass.
%	Length: L_u=5.291773e-11 m; Bohr radius
%	Energy: E_u=4.359827e-18 J; 27.2198eV
%	Time: t_u=2.41885e-17 s;
%
%	1 eV=1.6021917e-19;%1eV.
%	q_e=1.602191 C
%	eps0=8.854187e-12 [SI]
close all
clear all
clc
CC=sqrt(-1);%Complex imagnary unit.
%Parameters.
k=0;%Wave number. Try k=1,-1!
m=1;%Particle mass.
a=-10;%Beginning of interval.
b=10;%End of interval.
dx=0.5;%Standard deviation of x.
x0=(a+b)/2;%Average of x.
dt=1e-2;%Time step.
ntimes=200;%Number of time step.
render_dt=0.001; %Frame time interval.
N=64;%Number of sample points.
x=( [ 0 : 2*N-1 ] / (2*N) ) * (b-a) + a; %The sampling of x is very important

Psi0=exp(CC*k*x).*(2*pi*dx^2)^(-0.25).*exp(-(x-x0).^2/4/dx^2);%Initial wave function of the particle.
y_max=max(abs(Psi0).^2);

V=zeros(size(x)); %Potential function.

Psi_norm=trapz(x,Psi0.*conj(Psi0))
%Energy.
E=zeros(ntimes+1,1);
KO_Psi0= - nabla1d(Psi0,a, b,2)/2/m;
H_Psi0=KO_Psi0 + V.*Psi0;
E(1)=trapz(x,H_Psi0.*conj(Psi0));

%p
p=zeros(ntimes+1,1);
p_Psi0=-CC*nabla1d(Psi0,a, b,1);
p(1)=trapz(x,p_Psi0.*conj(Psi0));

%Time propagation
Psi=Psi0;

for i=1:ntimes
		t=i*dt;
		%Step 1 of SPO. exp(-i/h_bar/2 K dt)
		Psi=kinetic_energy_operator_exp(Psi,a, b, m, dt);
		
	  %Step 2 of SPO. exp(-i/h_bar/2 K dt)
   Psi=exp(-CC*V*dt) .* Psi;

   %Step 3 of SPO
   Psi=kinetic_energy_operator_exp(Psi,a, b, m, dt);
   %End nunerical solution for one time step (dt).

   if k==0
     %Analitical solution k = 0.
     Psi_anal=(2*pi*dx^2)^(-0.25).* ( 1+t.^2/(4*m^2*dx^4) ) .^(-0.25).* exp(-(x-x0).^2/4/dx^2/(1+t.^2/(4*m^2*dx^4)).*(1-CC.*t/2/m/dx^2));

     if i==1
       h_plot=plot(x,abs(Psi).^2, '+',x,abs(Psi_anal).^2);
       axis([-inf inf 0.01 y_max])
       xlabel( 'X')
       ylabel( 'PDF')
       legend( 'Numerical Solution' ,'Analytical Solution')
       drawnow
       pause(render_dt)
     else
       set(h_plot(1), 'YData',abs(Psi).^2)
       set(h_plot(2), 'YData',abs(Psi_anal).^2)
       pause(render_dt)
     end

   else

     if i==1
       h_plot2=plot(x,abs(Psi).^2);
       %axis([-inf inf 0.1 y_max])
       xlabel( 'X')
      ylabel( 'PDF')
      legend( 'Numerical')
      pause(render_dt)
	  else
		  set(h_plot2, 'YData',abs(Psi).^2);
		  pause(render_dt)
	  end


	end

	%Average of Energy
	KO_Psi=nabla1d(-Psi/2/m,a, b,2);
	H_Psi=KO_Psi + V.*Psi;
	E(i+1)=trapz(x,H_Psi.*conj(Psi));

	%Average of impulse. Operator of impusle: h_bar/i d/dx.
	p_Psi=-CC*nabla1d(Psi,a, b,1);
	p(i+1)=trapz(x,p_Psi.*conj(Psi));
end

figure(2)
subplot(211)
semilogy(1:ntimes+1,real(E))
xlabel('Time step')
ylabel('Average energy' )
subplot(212)
if k==0
  plot(1:ntimes+1,real(p))
else
  semilogy(1:ntimes+1,real(p))
end
xlabel('Time step')
ylabel('Average impulse' )