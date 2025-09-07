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

% measure difference to theoretical local bifurcation branch for square
% patterns with k = 2k_0 = 1, and thus M^*(2k_0) = 8.
Mstar = 8; % copy numerically found bifurcation point
k0 = 1;
c2 = -(g+k0^2)*(104 + 203*k0^2)/(24*pi^2*(1+k0^2+k0^4));

[grid,~,~] = getpte(p);
xi = transpose(2*(cos(k0*grid(1,:))+cos(k0*grid(2,:)))); 
quadcoeff = -2*(g+k0^2)*(104*g+203*k0^2)/(3*k0^2);
uapprox = sqrt((M-Mstar)/quadcoeff)*xi;
diffSqu = (p.mat.vM*((uSol-uapprox).^2))^(1/2);%p.mat.vM*(uSol - 2*uapprox).^2;

l2log = p.mat.vM*(log(1+uSol).^2);

out=[par;     % parameters 
    min(uSol);             % min u
    l2log;
    diffSqu];