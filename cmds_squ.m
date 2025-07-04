%% main command file for the bifurcation analysis of a 2D thermocapillary thin-film equation
% 
% There are two main command files:
% - cmds_hex is set up to detect hexagonal patterns
% - cmds_squ is set up to detect square patterns
%
% There are three additional command files to study patterns bifurcating at
% M=20.
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
p0=cswibra('init-squ','bpt1'); 
p=seltau(p0, 2,'squ',3); 
p.sol.ds=0.01;                      % -0.01 gives rise to "down-squares"
p=pmcont(p,206); 

%% c4: secondary bifurcation
p = swibra('squ','bpt3','square-sec-bif');  % bpt3 admissible pattern with 8 min around a max
p.sol.ds=0.01;
p=pmcont(p,121);

%% c5: plot bifurcation diagram squares
hold on;
plotbra('init-squ','cl','k');
plotbra('squ','cl','b');
plotbra('square-sec-bif','cl','b');
hold off;
if saveFigures
    set(gcf,'position',[0,0,680,400])
    saveas(gcf,'bifurcation-diag-squ','epsc');
end

%% c6: plot solution close to the bifurcation point
plotsol('squ','pt5',1,1,2,'cm','jet')   %third numeric argument: 1 mesh graph, 2 contour plot, 3 graph
if saveFigures
    set(gcf,'position',[0,0,400,400])
    saveas(gcf,'small-squ','epsc')
end

%% c7: plot solution close to film rupture
plotsol('squ','pt206',1,1,3,'cm','jet')
if saveFigures
    set(gcf,'position',[0,0,400,400])
    saveas(gcf,'film-rupture-squ','epsc')
end

%% c8: plot sec-bif solution
plotsol('squ','bpt3',1,1,2,'cm','jet')
if saveFigures
    set(gcf,'position',[0,0,400,400])
    saveas(gcf,'bpt3-squ','epsc')
end
plotsol('square-sec-bif','pt10',2,1,2,'cm','jet')
if saveFigures
    set(gcf,'position',[0,0,400,400])
    saveas(gcf,'sec-squ-pt10','epsc')
end
plotsol('square-sec-bif','pt25',3,1,2,'cm','jet')
if saveFigures
    set(gcf,'position',[0,0,400,400])
    saveas(gcf,'sec-squ-pt25','epsc')
end
plotsol('square-sec-bif','pt120',4,1,2,'cm','jet')
if saveFigures
    set(gcf,'position',[0,0,400,400])
    saveas(gcf,'film-rupture-sec-squ','epsc')
end

%% c9: plot integration constant \lambda = K(0) - K(v)
hold on;
plotbra('init-squ',15,3,'cl','k');
plotbra('squ',15,3,'cl','k');
plotbra('square-sec-bif',15,3,'cl','k');
hold off;
if saveFigures
    set(gcf,'position',[0,0,680,400])
    saveas(gcf,'constant-K-squ','epsc');
end

%% c10: plot L^2 norm of log(1+v)
hold on;
plotbra('init-squ',16,5,'cl','k');
plotbra('squ',16,5,'cl','k');
plotbra('square-sec-bif',16,5,'cl','k');
hold off;
if saveFigures
    set(gcf,'position',[0,0,680,400])
    saveas(gcf,'l2norm-squ','epsc');
end

%% c11: plot minimum
hold on;
plotbra('init-squ',17,4,'cl','k');
plotbra('squ',17,4,'cl','k');
plotbra('square-sec-bif',17,4,'cl','k');
hold off;
if saveFigures
    set(gcf,'position',[0,0,680,400])
    saveas(gcf,'minimum-squ','epsc');
end

%% c12: export of plot data
% extract and save solution data.

% load solution close to bifurcation point
data_squ_pt5 = loadp('squ','pt5');
sol_squ_pt5 = [data_squ_pt5.pdeo.grid.p;data_squ_pt5.u(1:end-3).'];

% load solution close to film-rupture
data_squ_pt206 = loadp('squ','pt206');
sol_squ_pt206 = [data_squ_pt206.pdeo.grid.p;data_squ_pt206.u(1:end-3).'];

% load secondary bifurcation solution
data_squ_bpt3 = loadp('squ','bpt3');
sol_squ_bpt3 = [data_squ_bpt3.pdeo.grid.p;data_squ_bpt3.u(1:end-3).'];

data_sec_squ_pt10 = loadp('square-sec-bif','pt10');
sol_sec_squ_pt10 = [data_sec_squ_pt10.pdeo.grid.p;data_sec_squ_pt10.u(1:end-3).'];

data_sec_squ_pt25 = loadp('square-sec-bif','pt25');
sol_sec_squ_pt25 = [data_sec_squ_pt25.pdeo.grid.p;data_sec_squ_pt25.u(1:end-3).'];

data_sec_squ_pt120 = loadp('square-sec-bif','pt120');
sol_sec_squ_pt120 = [data_sec_squ_pt120.pdeo.grid.p;data_sec_squ_pt120.u(1:end-3).'];

% extract bifurcation branch data
data_init = loadpp('init-squ');
brainit = data_init.branch([4,6],:);

data_squ = loadpp('squ');
bra_squ = data_squ.branch([4,6],:);

data_sec_squ = loadpp('square-sec-bif');
bra_sec_squ = data_sec_squ.branch([4,6],:);

% extract other plots data
data_init_const_squ = loadpp('init-squ');
bra_init_const_squ = data_init_const_squ.branch([4,9],:);

data_init_l2norm_squ = loadpp('init-squ');
bra_init_l2norm_squ = data_init_l2norm_squ.branch([4,11],:);

data_init_min_squ = loadpp('init-squ');
bra_init_min_squ = data_init_min_squ.branch([4,10],:);


data_const_squ = loadpp('squ');
bra_const_squ = data_const_squ.branch([4,9],:);

data_l2norm_squ = loadpp('squ');
bra_l2norm_squ = data_l2norm_squ.branch([4,11],:);

data_min_squ = loadpp('squ');
bra_min_squ = data_min_squ.branch([4,10],:);


data_sec_const_squ = loadpp('square-sec-bif');
bra_sec_const_squ = data_sec_const_squ.branch([4,9],:);

data_init_l2norm_squ = loadpp('square-sec-bif');
bra_sec_l2norm_squ = data_init_l2norm_squ.branch([4,11],:);

data_init_min_squ = loadpp('square-sec-bif');
bra_sec_min_squ = data_init_min_squ.branch([4,10],:);

% write data to text files
if saveData
    formatSpec = '(%5.5f,%5.5f) ';

    fileID = fopen('branchDataSqu.txt','w');
    fprintf(fileID,'%1s\n\n',['==================',' brainit ','==================']);
    fprintf(fileID,formatSpec,brainit);

    fprintf(fileID,'\n\n%1s\n\n',['==================',' bra_squ ','==================']);
    fprintf(fileID,formatSpec,bra_squ);

    fprintf(fileID,'\n\n%1s\n\n',['==================',' bra_sec_squ ','==================']);
    fprintf(fileID,formatSpec,bra_sec_squ);
    fclose(fileID);


    fileID = fopen('solutionDataSqu.txt','w');
    fprintf(fileID,'%1s\n\n',['==================',' solution squares pt5 ','==================']);
    fprintf(fileID,formatSpec,sol_squ_pt5);

    fprintf(fileID,'\n\n%1s\n\n',['==================',' solution squares pt206 ','==================']);
    fprintf(fileID,formatSpec,sol_squ_pt206);
    fclose(fileID);

    
    fileID = fopen('solutionDataSquSec.txt','w');
    fprintf(fileID,'%1s\n\n',['==================',' solution squares secondary bif pt ','==================']);
    fprintf(fileID,formatSpec,sol_squ_bpt3);

    fprintf(fileID,'\n\n%1s\n\n',['==================',' solution squares secondary branch pt10 ','==================']);
    fprintf(fileID,formatSpec,sol_sec_squ_pt10);

    fprintf(fileID,'\n\n%1s\n\n',['==================',' solution squares secondary branch pt25 ','==================']);
    fprintf(fileID,formatSpec,sol_sec_squ_pt25);

    fprintf(fileID,'\n\n%1s\n\n',['==================',' solution squares secondary branch pt120 ','==================']);
    fprintf(fileID,formatSpec,sol_sec_squ_pt120);
    fclose(fileID);

    
    fileID = fopen('constantPlotSqu.txt','w');
    fprintf(fileID,'%1s\n\n',['==================',' constant_init ','==================']);
    fprintf(fileID,formatSpec,bra_init_const_squ);

    fprintf(fileID,'\n\n%1s\n\n',['==================',' constant ','==================']);
    fprintf(fileID,formatSpec,bra_const_squ);

    fprintf(fileID,'\n\n%1s\n\n',['==================',' constant_sec ','==================']);
    fprintf(fileID,formatSpec,bra_sec_const_squ);
    fclose(fileID);
    
    
    fileID = fopen('l2normPlotSqu.txt','w');
    fprintf(fileID,'%1s\n\n',['==================',' l2norm_init ','==================']);
    fprintf(fileID,formatSpec,bra_init_l2norm_squ);

    fprintf(fileID,'\n\n%1s\n\n',['==================',' l2norm ','==================']);
    fprintf(fileID,formatSpec,bra_l2norm_squ);

    fprintf(fileID,'\n\n%1s\n\n',['==================',' l2norm_sec ','==================']);
    fprintf(fileID,formatSpec,bra_sec_l2norm_squ);
    fclose(fileID);


    fileID = fopen('minSqu.txt','w');
    fprintf(fileID,'%1s\n\n',['==================',' min_init ','==================']);
    fprintf(fileID,formatSpec,bra_init_min_squ);

    fprintf(fileID,'\n\n%1s\n\n',['==================',' min ','==================']);
    fprintf(fileID,formatSpec,bra_min_squ);

    fprintf(fileID,'\n\n%1s\n\n',['==================',' min_sec ','==================']);
    fprintf(fileID,formatSpec,bra_sec_min_squ);
    fclose(fileID);
end