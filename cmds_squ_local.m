%% main command file for the bifurcation analysis of a 2D thermocapillary thin-film equation
% 
% This is an additional command file to study compare the local behaviour
% of the numerical continuation branch of square patterns with the rigorous
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
lx = [2*pi 2*pi];                   % set domain size; this ratio does not allow for hexagonal pattern
nx = [100 100];                     % number of discretisation points per dimension in domain
Minit = 7.9;                        % initial Marangoni number
ginit = 1;                          % set gravitational constant
lambdaInit = 0;                     % set initial integration constant (lambda = M*( K(0) - K ))
par = [Minit, ginit, lambdaInit];
p = tfinit(p,lx,nx,par);
p.fuha.outfu = @tfbra;
p = setfn(p,'init-squ');
para = 1;                           % set Marangoni number as bifurcation parameter
p.nc.dsmax=0.03;
                                                                                                    
%% c2: continuation of the trivial branch
p = cont(p,15);                     % set up to detect bifurcation point at M = 8 (two full squares)

%% c3: switch branch to periodic bifurcation branches and continue up to film-rupture through squares
p.nc.dsmax=0.0001;
p0=cswibra('init-squ','bpt1'); 
p=seltau(p0, 2,'squ',3); 
p.sol.ds=0.0001;                     
p=pmcont(p,25); 

%% c4: plot difference of numerical solution and local approximation
figure(16);
hold on;
plotbra('squ',16,6,'cl','k');
hold off;
if saveFigures
    set(gcf,'position',[0,0,680,400])
    saveas(gcf,'error-squ-max','epsc');
end

%% c5: plot approximation and numerical branch
figure(17);
hold on;
plotbra('squ',17,8,'cl','k'); % theoretical prediction in black
plotbra('squ','wnr',17,'cl','b'); % numerical branch in blue
hold off;
if saveFigures
    set(gcf,'position',[0,0,680,400])
    saveas(gcf,'comparisonSqu','epsc');
end

%% c12: export of plot data
% extract and save solution data.

data_error_squ = loadpp('squ');
bra_error_squ = data_error_squ.branch([4,12],:);

data_squ = loadpp('squ');
bra_squ = data_squ.branch([4,6],:);

data_approx_branch_squ = loadpp('squ');
bra_approx_branch_squ = data_approx_branch_squ.branch([4,14],:);

% write data to text files
if saveData
    formatSpec = '(%5.5f,%5.5f) ';
    
    fileID = fopen('comparisonSqu.txt','w');
    fprintf(fileID,'\n\n%1s\n\n',['==================',' bra_squ ','==================']);
    fprintf(fileID,formatSpec,bra_squ);

    fprintf(fileID,'\n\n%1s\n\n',['==================',' bra_approx_squ ','==================']);
    fprintf(fileID,formatSpec,bra_approx_branch_squ);
    fclose(fileID);

    fileID = fopen('errorSqu.txt','w');
    fprintf(fileID,'\n\n%1s\n\n',['==================',' bra_error_squ ','==================']);
    fprintf(fileID,formatSpec,bra_error_squ);
    fclose(fileID);
end