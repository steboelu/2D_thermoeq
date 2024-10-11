function r = sG(p,u)
% SG evaluates bifurcation function
% F(v) = \Lambda v - f(v)
% with f(v) = gv - M(log((1+v)/(2+v))+1/(2+v))+M*(K(0)+lambda)
%%%%%%%%%%%%%%%%%%%% TYPO notation for constants not consistent
%
%   r = sG(p,u)
%
% Input:
%   p : 'p' variable of pde2path (contains collection of current
%       continuation data
%   u : solution (entries 1:p.np) and parameters (entries (p.nu+1:end))
%
% Ouput:
%   r : 'remainder' of equation F(v) = 0
%

f = nodalf(p,u); % evaluate f(v)
K = p.mat.K; % stiffness matrix from FEM descritization
r = K*u(1:p.nu) - p.mat.M*f; % calculate residual
end