function Gu = sGjac(p,u)
% SGJAC evaluates v-derivative of bifurcation function
% 
%   F(v) = \Lambda v - f(v)     with
%   f(v) = gv - M*( 1/(2+v) + log((1+v)/(2+v)) ) + M*K(0) - lambda
%
%   Gu = sGjac(p,u)
%
% Input:
%   p : 'p' variable of pde2path (contains collection of current
%       continuation data
%   u : solution (entries 1:p.np) and parameters (entries (p.nu+1:end))
%
% Output:
%   Gu : derivative of F(v)
%

par = u(p.nu+1:end);                    % par = [M,g,lambda]
M = par(1);                             % Marangoni number
g = par(2);                             % gravity
n = p.np;                               % number of grid points
u1 = u(1:n);                            % extract solution u

fu = (-g + M*(1./((1+u1).*(2+u1).^2))); % calculate v-derivative of 'reaction term' f(v)
Fu = spdiags(fu,0,n,n);                 % set up sparse matrix with diagonal entries f(u)
K = p.mat.K;                            % stiffness matrix of FEM discretization
Gu = K-p.mat.M*Fu;                      % calculate derivative
end