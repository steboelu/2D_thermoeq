function q=qf(p,u)
% QF help functions to calculate the current mass (integral over one
% period) for the mass constraint
%
%   q = qf(p,u)
%
% Input:
%   p : 'p' variable of pde2path (contains collection of current
%       continuation data
%   u : solution (entries 1:p.np) and parameters (entries (p.nu+1:end))
%
% Output:
%   q : current mass
%

u=u(1:p.nu);                    % extract solution
q=p.mat.vM*u/p.vol;             % calculate mass
end