function [A,r] = laplacian(order,n,h,test)
% [A,r] = laplacian(n,h,order,test)
%
% Input:
% order : Interior order of accuracy
% n     : Number of grid points 
% h     : Grid spacing
% 
% Optional:
%  test (true/false) : Tests that the SBP operators are accurate.
%
% Output:
% A     : Spatial discretization (n x n Sparse matrix)
% rm    : Radial coordinates (n x 1 vector)
%
% Discretization of the diffusion term in radial coordinates
% v_t = 1/r*(r*v_r)_r
% Spatial discretization:
% v_t = A*v,

if nargin < 4
  test = false;
end

truncate = true;
[rp,rm,Pp,Pm,Qp,Qm] = sbp_staggered(order,n,h,truncate,test);
r  = rm;

% Assemble radial coordinate matrices
Rp = spdiags(rp,0,n+1,n+1);
Rm = spdiags(rm,0,n,n);

% Assemble Laplacian
A = inv(Rm)*inv(Pm)*Qm*Rp*inv(Pp)*Qp;
