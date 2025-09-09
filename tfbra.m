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
% Mstar = 8.00090265751118; % grid 200x200
Mstar = 8.00359895163406; % grid 100x100
k0 = 1;
[grid,~,~] = getpte(p);

c2 = -(g+k0^2)*(104 + 203*k0^2)/(24*pi^2*(1+k0^2+k0^4));
hatxi = transpose(cos(k0*grid(1,:))+cos(k0*grid(2,:)));
normxi = 8*pi/k0;
hatxi = hatxi/(p.mat.vM*hatxi.^2)^(1/2);
uapprox1 = -sqrt((M-Mstar)/c2)*hatxi;

xi = transpose(2*(cos(k0*grid(1,:))+cos(k0*grid(2,:)))); 
quadcoeff = -2*(g+k0^2)*(104*g+203*k0^2)/(3*k0^2);
uapprox = -sqrt((M-Mstar)/quadcoeff)*xi;
diffSqu = max(abs(uSol-uapprox));

diffSqu = (p.mat.vM*(uapprox.^2))^(1/2);

diffSqu2 = (p.mat.vM*(uapprox1.^2))^(1/2);

diffSqu3 = (p.mat.vM*(sqrt((M-Mstar)/c2)*(xi/normxi)).^2)^(1/2);



% ii = 1;
% if (ii < 2 & M < 8)
%     figure(30);
%     [X,Y] = meshgrid(grid(1,:),grid(2,:));
%     xi = 2*transpose(cos(k0*X)+cos(k0*Y));
%     uapprox = sqrt((M-Mstar)/quadcoeff)*xi;
%     surfc(X,Y,uapprox);
%     axis off;
%     grid off;
%     view(2);
%     drawnow;
%     pause();
%     ii = 2;
% end

l2log = p.mat.vM*(log(1+uSol).^2);

out=[par;     % parameters 
    min(uSol);             % min u
    l2log;
    diffSqu;
    diffSqu2;
    diffSqu3];