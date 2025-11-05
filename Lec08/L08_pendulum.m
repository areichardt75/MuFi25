%% Mufi 2025 - Lecture 8. - Solution of ODEs
% Simple pendulum 

% general form : Th'' = -g/ell*sin(Th)
% linearized form : Th'' = -g/ell*Th

ell = 1;
g = 9.81;

% odefun : x = [Th Th']
odefun = @(x,ell,g) [x(2);-g/ell*sin(x(1))];
odelin = @(x,ell,g) [x(2);-g/ell*x(1)];

T0 = 2*pi/sqrt(g/ell);
Tspan = [0 2.5*T0];
x0 = [pi/3;0];
[T,Y] = ode45(@(t,x) odefun(x,ell,g), Tspan, x0);
[Tlin,Ylin] = ode45(@(t,x) odelin(x,ell,g), Tspan, x0);

Th = Y(:,1); Thv = Y(:,2);
Thlin=Ylin(:,1); Thvlin=Y(:,2);

%% Visualization
figure; 
  plot(T,Th,'r-',Tlin,Thlin,'k-','LineWidth',2);
  xlabel('time'); 
  ylabel('Theta');
  legend('general','linearized','Location','northeast');