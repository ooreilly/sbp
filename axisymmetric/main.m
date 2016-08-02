% --------------------------------------------------------- 
% Solves the problem:
% v_t = 1/r*(r*v_r)_r, 
% with a dirichlet boundary condition is specified at r = R.
% The numerical solution is determined in the interval
% h/2 <= r <= R - h/2, 
% where h is the grid spacing. 
% --------------------------------------------------------- 

% order: Order of accuracy (2,4,6, or 8)
%     R: Domain size
%     n: Number of grid points
%     h: Grid spacing
%    v0: Initial condition

order = 8;
R     = 20;
n     = 60;
h     = R/n;
v0    = @(r) exp(-(r-R/2).^2);

[A r] = laplacian(order,n,h);
f = @(t,v) A*v;

[T,V] = ode45(f,[0 10],v0(r));

for i=1:length(T)
  plot(r,V(i,:),'bo-')
  xlabel('r');
  ylabel('v');
  ylim([0 1]);
  drawnow;
end

