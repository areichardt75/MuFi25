% MuFi25 / Lecture 10. Hyperbolic type PDE
% -------------------------------------------
% Solution of wave-equation using a Lax-Wendroff one-step method
% 
% f,g 
% fi(x) = f(x,0) = g(x,0)
% c = u*dt/dx

L = 20;
dx = 0.1;
x = 0:dx:L;
xlen = length(x);

Tmax = 8;
dt = 0.1;
t = 0:dt:Tmax;
tlen = length(t);

u = 0.5; 
c = u*dt/dx; 
fprintf('c = %8.3f \n',c);

f = zeros(tlen, xlen);
g = zeros(tlen, xlen);
f(1,:) = phi(x);
% g(1,:) = f(1,:)/2;

for id=2:tlen
  fprev = f(id-1,:);
  gprev = g(id-1,:);
  [fnext, gnext] = LWone(fprev, gprev, c);
  f(id,:) = fnext;
  g(id,:) = gnext;
end

%%
tidx = [1 floor([tlen/5 tlen/3 tlen/2 3*tlen/4 4*tlen/5]) tlen];
plotu(x,f,tidx);


%% 
plotu(x,f,[1 floor(tlen/3) floor(tlen/2) floor(tlen/3*2) tlen])

%% Functions 

function [fnext, gnext] = LWone(fprev, gprev, c)
  fnext = zeros(size(fprev));
  gnext = zeros(size(gprev));
  for id=2:length(fprev)-1
    fnext(id) = fprev(id)-0.5*c*(gprev(id+1)-gprev(id-1))+...
      c*c/2*(fprev(id+1)-2*fprev(id)+fprev(id-1));
    gnext(id) = gprev(id)-0.5*c*(fprev(id+1)-fprev(id-1))+...
      c*c/2*(gprev(id+1)-2*gprev(id)+gprev(id-1));
  end
  fnext(1) = fprev(1); fnext(end)=fprev(end);
  gnext(1) = gprev(1); gnext(end) = gprev(end);

end



%% - ------------
% Other (utility) functions 
function figh = plotu(x,u,tpoints)
  figh = figure; 
  id = tpoints(1);
  plot(x,u(id,:),'r-','LineWidth',2); xlabel('x'); ylabel('T(x)');
  fprintf('plotting : %3d --> %3d | ',id,tpoints(id));

  if length(tpoints)>1
    hold on;
    for id=2 : length(tpoints)
      plot(x,u(tpoints(id),:),'-','LineWidth',2);
      fprintf('%3d --> %3d | ',id,tpoints(id));
    end
  end
  fprintf('\n');
end

function uIC = phi(x)
%   uIC = zeros(size(x));
  x0 = 5; 
  d = 0.5;
  umax = 100;
%   uIC = 200*((x>=5).*(x-5).*(x<=5.5)+(x>5.5).*(1-(x-5)).*(x<=6));
  uIC = 100*(x>=5).*sin(2*pi*(x-5)/5).*(x<=10);
%   uIC = 2*umax*((x>=(x0-d)).*(x-x0-d).*(x<=x0)+...
%     (x>x0).*(1-(x-x0)).*(x<=(x0+d)));
end


