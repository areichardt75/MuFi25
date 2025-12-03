
% Lax Method for solution of 1D convection problem
% 

% Setting parameters
L = 20;
dx = 0.2;   % equidistant grid
x = 0:dx:L; % creating grid
Nx = length(x);

% 
uv = 1;
dt = 0.002;
Tmax = 20;
c = uv*dt/dx;

t = 0:dt:Tmax;
u = zeros(length(t),length(x));
u(1,:) = setIC(x);
uLW = u;
uBTCS = u;
% plotu(x,u,1);
%% 
for id=2:length(t)
%   un = Lax(u(id-1,:),c);
%   u(id,:) = un;
  unLW = LWonestep( uLW(id-1,:), c);
  uLW(id,:)=unLW;
%   fprintf('/');
  unBTCS = BTCS(uBTCS(id-1,:),c);
%   fprintf('\\');
  uBTCS(id,:) = unBTCS;
%   fprintf('-');
end

%% 
tidx = [1 10 20 length(t)];
% f1 = plotu(x,u,tidx);
%% 
f1LW=plotu(x,uLW,tidx);
title('Lax-Weindroff');
legend('1','10','20','end','location','best');
%% Backward Time scheme
[f1BT,hpBT] = plotu(x, uBTCS, tidx);
% f2 = plotu(x,u,40);
title('Backward Time Centered Space method');
legend('1','10','20','end','location','best');

%% 3D plot / 
hh = plot3du(x,t,uLW);

%% Functions 

function unext = Lax(uprev, c)
  unext = zeros(size(uprev));
  for id=2:length(uprev)-1
    unext(id) = 0.5*(uprev(id-1)+uprev(id+1))-0.5*c*(uprev(id+1)-uprev(id-1));
  end
  unext(1)=uprev(1);
  unext(end) = uprev(end);
end

function unext = LWonestep(uprev,c)
  unext = zeros(size(uprev));
  for id=2:length(uprev)-1
    unext(id) = uprev(id)-0.5*c*(uprev(id+1)-uprev(id-1))+0.5*c^2*(uprev(id+1)-2*uprev(id)+uprev(id-1));
  end
  unext(1) = uprev(1);
  unext(end) = uprev(end);

end


function unext = BTCS(uprev, c)
  A = zeros(length(uprev),length(uprev));
  B = transpose(uprev);
%   fprintf('|');
  for id=2:length(uprev)-1
    A(id,id-1:id+1) = [-c/2 1 c/2];
  end
  A(1,1) = 1;
  A(end,end)=1;
  unext = A \ B;
  unext = transpose(unext);
end

% Other (utility) functions 
function [figh,hplot] = plotu(x,u,tpoints)
  figh = figure; 
  id = tpoints(1);
  hplot = plot(x,u(id,:),'r-o'); xlabel('x'); ylabel('T(x)');
  fprintf('plotting : %3d --> %3d | ',id,tpoints(id));

  if length(tpoints)>1
    hold on;
    for id=2 : length(tpoints)
      plot(x,u(tpoints(id),:),'-d');
      fprintf('%3d --> %3d | ',id,tpoints(id));
    end
  end
  fprintf('\n');
end

function uIC = setIC(x)
%   uIC = zeros(size(x));
  x0 = 5; 
  d = 0.5;
  umax = 100;
%   uIC = 200*((x>=5).*(x-5).*(x<=5.5)+(x>5.5).*(1-(x-5)).*(x<=6));
  uIC = 100*(x>=5).*sin(2*pi*(x-5)/5).*(x<=10);
%   uIC = 2*umax*((x>=(x0-d)).*(x-x0-d).*(x<=x0)+...
%     (x>x0).*(1-(x-x0)).*(x<=(x0+d)));
end

%% plot3du(x, t, u)
function hh = plot3du(x, t, u)
% function plot3du(x, t, u)
  [xx,tt] = meshgrid(x,t);
  figure;
    hh= surf(xx,tt, u);
    xlabel('x');
    ylabel('t'); 
    zlabel('u');
    view(30,65);
    hh.LineStyle = 'none';
end