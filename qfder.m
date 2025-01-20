function qu=qfder(p,u)
% QFDER help function to calculate d_u qf(u), where qf(u) denotes the mass
% of u
%
%   qu = qfder(p,u)
%
% Input:
%   p : 'p' variable of pde2path (contains collection of current
%       continuation data
%   u : solution (entries 1:p.np) and parameters (entries (p.nu+1:end))
%
% Output:
%   qu : derivative of mass function
%

qu=(1/p.vol)*p.mat.vM;
end