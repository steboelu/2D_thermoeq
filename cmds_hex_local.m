%% main command file for the bifurcation analysis of a 2D thermocapillary thin-film equation
% 
% This is an additional command file to study compare the local behaviour
% of the numerical continuation branch of up-hexagons with the rigorous
% prediction from local bifurcation theory.
%
% The command file performs a numerical continuation of the second order
% equation
%
%   0 = \Delta v - gv + M*( 1/(2+v) + log((1+v)/(2+v)) ) - M*K(0) + lambda
%
% Here, 1+v is the fluid height, g is the gravitation constant, M is the
% Marangoni number (the primary bifurcation parameter), and lambda is a
% constant to conserve the mass of the bifurcating solution.
%
% The numerical continuation is performed with pde2path with Neuman
% boundary conditions on a domain (-lx,lx)x(-ly,ly). 

%% clear mainspace (keep paths for pde2path) and close all plot windows
clc;
close all;
keep pphome;

% set to 1 to save plotted figures as eps
saveFigures = 1;
% set to 1 to save plot data
saveData = 1;

%% c1: init and set parameters
p = [];
lx = [2*pi (2/sqrt(3))*pi];         % set domain size; this ratio allows for hexagonal pattern
nx = [100 100];                     % number of discretisation points per dimension in domain
Minit = 7.8;                        % initial Marangoni number
ginit = 1;                          % set gravitational constant
lambdaInit = 0;                     % set initial integration constant (lambda = M*( K(0) - K ))
par = [Minit, ginit, lambdaInit];   
p = tfinit(p,lx,nx,par);
p.fuha.outfu = @tfbra;
p = setfn(p,'init-hex');
para = 1;                           % set Marangoni number as bifurcation parameter
p.nc.dsmax=0.03;
                                                                                                    
%% c2: continuation of the trivial branch
p = cont(p,15);                     % set up to detect bifurcation point at M = 8 (one full hexagon)

 %% c3: switch branch to periodic bifurcation branches and continue up to film-rupture through up-hexagons
p0=qswibra('init-hex','bpt1');
p=gentau(p0,[0.1,0.1],'hex-up');        % Detected tangent directions are \phi_1 = cos(k_1*(x,y)) and \phi_2 = cos(k_2*(x,y))+cos(k_3*(x,y)). 
%                                     Generate tanget direction \phi_1 + \phi_2, which yields hexagons.
p.sol.ds=0.001;                     % positive for up-hexagons
p=pmcont(p,25);

%% c4: plot the error between numerical continuation and rigorous prediction
plotbra('hex-up',30,7,'cl','k');
if saveFigures
    set(gcf,'position',[0,0,680,400])
    saveas(gcf,'error-hex-up-max','epsc');
end

%% c5: plot numerical branch and theoretical prediction
figure(31);
hold on;
plotbra('hex-up',31,9,'cl','k'); % theoretical prediction in black
plotbra('hex-up','wnr',31,'cl','b'); % numerical branch in blue
hold off;

if saveFigures
    set(gcf,'position',[0,0,680,400])
    saveas(gcf,'comparisonHexUp','epsc');
end

%% c6: export of plot data
% extract and save solution data.

data_error_hex_up = loadpp('hex-up');
bra_error_hex_up = data_error_hex_up.branch([4,13],:);

data_hex_up = loadpp('hex-up');
bra_hex_up = data_hex_up.branch([4,6],:);

data_approx_branch_hex_up = loadpp('hex-up');
bra_approx_branch_hex_up = data_approx_branch_hex_up.branch([4,15],:);

% write data to text files
if saveData
    formatSpec = '(%5.5f,%5.5f) ';
    
    fileID = fopen('comparisonHexUp.txt','w');
    fprintf(fileID,'\n\n%1s\n\n',['==================',' bra_hex_up ','==================']);
    fprintf(fileID,formatSpec,bra_hex_up);

    fprintf(fileID,'\n\n%1s\n\n',['==================',' bra_approx_hex_up ','==================']);
    fprintf(fileID,formatSpec,bra_approx_branch_hex_up);
    fclose(fileID);

    fileID = fopen('errorHexUp.txt','w');
    fprintf(fileID,'\n\n%1s\n\n',['==================',' bra_error_hex_up ','==================']);
    fprintf(fileID,formatSpec,bra_error_hex_up);
    fclose(fileID);
end