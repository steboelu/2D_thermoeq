function f = nodalf(p,u)
% NODALf help function, which calculates the 'reaction' part of the ode
% (see output)
%
%   f = nodal(p,u)
%
% Input:
%   p : 'p' variable of pde2path (contains collection of current
%       continuation data
%   u : solution (entries 1:p.np) + parameters (entries (p.nu+1:end))
%
% Output:
%   f : evaluation of -f(v) = -gv + M*( 1/(2+v) + log((1+v)/(2+v)) ) - M*K(0) + lambda
%

uSol = u(1:p.np);                       % extract solution and parameter data
par = u(p.nu+1:end);                    % par = [M,g,lambda]
g = par(2);                             % gravitational constant
M = par(1);                             % Marangoni number
lambda = par(3);                        % mass constraint parameter

K = -M*(log(1/2)+1/2)*ones(size(uSol)); % this is actually -M*K(0) in the notation of the paper.
f = -g*uSol+M*(log((1+uSol)./(2+uSol))+1./(2+uSol)) + K + lambda;   % this is actually -f(v) in the notation of the paper.
% recall lambda = M*( K(0) - K(v) ). Thus K + lambda = - MK(v) in the notation of the paper.
end