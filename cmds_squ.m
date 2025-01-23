%% main command file for the bifurcation analysis of a 2D thermocapillary thin-film equation
% 
% There are two main command files:
% - cmds_hex is set up to detect hexagonal patterns
% - cmds_squ is set up to detect square patterns
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
p = setfn(p,'init-squ');
para = 1;                           % set Marangoni number as bifurcation parameter
p.nc.dsmax=0.03;
                                                                                                    
%% c2: continuation of the trivial branch
p = cont(p,15);                     % set up to detect bifurcation point at M = 8 (two full squares)

%% c3: switch branch to periodic bifurcation branches and continue up to film-rupture through squares
p0=cswibra('init-squ','bpt1'); 
p=seltau(p0, 2,'squ',3); 
p.sol.ds=0.01; 
p=pmcont(p,300); 

%% c-sec: secondary bifurcation
p = swibra('squ','bpt3','square-sec-bif');
p.sol.ds=0.01;
p=pmcont(p,100);

% plot secondary bifurcations
if saveFigures
    set(gcf,'position',[0,0,500,400])
    plotsol('square-sec-bif','pt100',1,1,3,'cm','jet')
end
saveas(gcf,'squ-secondary-bif','epsc')

%% c4: plot bifurcation diagram squares
hold on;
plotbra('init-squ','cl','k');
plotbra('squ','cl','b');
hold off;
% save bifurcation diagram as eps
if saveFigures
    set(gcf,'position',[0,0,500,400])
    saveas(gcf,'bifurcation-diag-squ','epsc');
end

%% c5: plot solution close to the bifurcation point
if saveFigures
    set(gcf,'position',[0,0,500,400])
    plotsol('squ','pt5',1,1,2,'cm','jet')
end
saveas(gcf,'small-squ','epsc')

%% c6: plot solution close to film-rupture
if saveFigures
    set(gcf,'position',[0,0,500,400])
    plotsol('squ','pt120',1,1,2,'cm','jet')
end
saveas(gcf,'film-rupture-squ','epsc')

%% c7: export of plot data
% extract and save solution data.

% load solution close to bifurcation point
data_squ_pt5 = loadp('squ','pt5');
sol_squ_pt5 = [data_squ_pt5.pdeo.grid.p;data_squ_pt5.u(1:end-3).'];

% load solution close to film-rupture
data_squ_pt202 = loadp('squ','pt202');
sol_squ_pt202 = [data_squ_pt202.pdeo.grid.p;data_squ_pt202.u(1:end-3).'];

% extract bifurcation branch data
data_init = loadpp('init-squ');
brainit = data_init.branch([4,6],:);

data_squ = loadpp('squ');
bra_squ = data_squ.branch([4,6],:);

% write data to text files
if saveData
    formatSpec = '(%5.5f,%5.5f) ';

    fileID = fopen('branchData.txt','w');
    fprintf(fileID,'%1s\n\n',['==================',' brainit ','==================']);
    fprintf(fileID,formatSpec,brainit);

    fprintf(fileID,'\n\n%1s\n\n',['==================',' bra_squ ','==================']);
    fprintf(fileID,formatSpec,bra_squ);
    fclose(fileID);

    fileID = fopen('solutionData.txt','w');
    fprintf(fileID,'%1s\n\n',['==================',' solution squares pt5 ','==================']);
    fprintf(fileID,formatSpec,sol_squ_pt5);
    fclose(fileID);

    fileID = fopen('solutionData.txt','w');
    fprintf(fileID,'%1s\n\n',['==================',' solution squares pt202 ','==================']);
    fprintf(fileID,formatSpec,sol_squ_pt202);
    fclose(fileID);
end