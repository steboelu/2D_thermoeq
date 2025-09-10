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

%% measure difference to theoretical local bifurcation branch for square

[grid,~,~] = getpte(p); % get grid

% patterns with k = 2k_0 = 1, and thus M^*(2k_0) = 8.
% Mstar = 8.00090265751118; % grid 200x200
Mstar = 8.00359895163406; % grid 100x100
% Mstar = 5.00004327981216; % grid 200x200 for k0 = 1/2
k0 = 1;

c2 = -(g+k0^2)*(104 + 203*k0^2)/(24*pi^2*(1+k0^2+k0^4)); % M''(0) for squares
hatxi = 2*transpose(cos(k0*grid(1,:))+cos(k0*grid(2,:))); % non-normalised kernel element
normxi = 4*pi*(k0^2+1+1/k0^2)^(1/2); % norm of kernel element
hatxi = hatxi/normxi; % normalised kernel element
uapprox1 = sqrt(2*(M-Mstar)/c2)*hatxi; % theoretical prediction

SquApproxL2 = (p.mat.vM*(uapprox1.^2))^(1/2); % L^2-norm of theoretical prediction
diffSqu = max(abs(uSol-uapprox1)); % L^\infty-norm of difference between numerics and theory


% measure difference between local bifucation in hexagonal case
% patterns with k = k0 = 1 and thus, M^*(k_0) = 8
Mstar = 8.00110018265581; % copy numerically found bifurcation point Mstar, grid 100x100
k0 = 1;
hatpsi = 2*transpose(cos(k0*grid(1,:)) + cos(k0*(-grid(1,:)/2+sqrt(3)*grid(2,:)/2)) + cos(k0*(-grid(1,:)/2-sqrt(3)*grid(2,:)/2))); % non-normalised kernel element
normpsi = 4*3^(1/4)*sqrt((1+k0^2+k0^4)/k0^2)*pi; % norm of kernel element
alpha = 2*(g+k0^2)/(3^(1/4)*pi*sqrt(k0^2+1+1/k0^2)); % M'(0) for hexagons
uapproxHex = ((M - Mstar)/alpha)*hatpsi/normpsi; % theoretical prediction

HexApproxL2 = (p.mat.vM*(uapproxHex.^2))^(1/2); % L^2 norm of theoretical prediction

diffHex = max(abs(uapproxHex-uSol)); % L^\infty-norm of difference between numerics and theory

% figure(30)
% p.pdeo.grid.plot(uapproxHex-uSol,'LineStyle','none'); view(2);
% colorbar;
% pause();


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
    diffHex;
    SquApproxL2;
    HexApproxL2];