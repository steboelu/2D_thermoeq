function out=tfbra(p,u)
% output to bifurcation diagram function for fch 
% parameters, min u, L^2-norm of log(u+1)


uSol=u(1:p.nu);                                                         % extract solution
utriv = zeros(size(uSol));                                              % trivial solution
par = u(p.nu+1:end);                                                    % par = [M,g,lambda]
g = par(2);                                                             % gravity
M = par(1);                                                             % Marangoni number
% E = p.mat.vM*((1/2)*(1/wnr)^2*(p.mat.Dx*uSol).^2+g*(1+uSol).^2 - M*(1+uSol).*log((1+uSol)./(2+uSol))); % calculate the energy
% deltaE = E - p.mat.vM*(g*(1+utriv).^2 - M*(1+utriv).*log((1+utriv)./(2+utriv)));

l2log = p.mat.vM*(log(1+uSol).^2);

out=[par;     % parameters 
    min(uSol);             % min u
    l2log];