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
nx = [200 200];                     % number of discretisation points per dimension in domain
Minit = 19.9;                       % initial Marangoni number
ginit = 1;                          % set gravitational constant
lambdaInit = 0;                     % set initial integration constant (lambda = M*( K(0) - K ))
par = [Minit, ginit, lambdaInit];
p = tfinit(p,lx,nx,par);
p.fuha.outfu = @tfbra;
p = setfn(p,'init-squ20');
para = 1;                           % set Marangoni number as bifurcation parameter
p.nc.dsmax=0.03;
                                                                                                    
%% c2: continuation of the trivial branch
p = cont(p,15);                     

%% c3: switch branch to periodic bifurcation branches and continue up to film-rupture through squares
p0=cswibra('init-squ20','bpt1'); 
p=seltau(p0, 1,'squ20',3); 
p.sol.ds=0.01;                     
p=pmcont(p,500); 

%% c4: secondary bifurcation
% bpt3, bpt6 admissible patterns
p = swibra('squ20','bpt3','square20-sec-bif3');
p.sol.ds=0.01;
p=pmcont(p,400);
p = swibra('squ20','bpt6','square20-sec-bif6');
p.sol.ds=0.01;
p=pmcont(p,400);

%% c5: plot bifurcation diagram squares M=20
hold on;
plotbra('init-squ20','cl','k');
plotbra('squ20','cl','b');
plotbra('square20-sec-bif6','cl','b');
plotbra('square20-sec-bif3','cl','b');
hold off;
if saveFigures
    set(gcf,'position',[0,0,680,400])
    saveas(gcf,'bifurcation-diag-squ20','epsc');
end

%% c6: plot joined diagram: primary branches M=8, M=20
hold on;
plotbra('init-squ20','cl','k');
plotbra('squ20','cl','b');
plotbra('init-squ','cl','k');
plotbra('squ','cl','b');
hold off;
if saveFigures
    set(gcf,'position',[0,0,680,400])
    saveas(gcf,'bifurcation-diag-squ-joined','epsc');
end

%% c7: plot solution close to film rupture
plotsol('squ20','pt500',1,1,2,'cm','jet');
if saveFigures
    set(gcf,'position',[0,0,400,400])
    saveas(gcf,'large-squ20','epsc')
end

%% c8: plot solution from secondary bpt6 close to film rupture
plotsol('square20-sec-bif6','pt198',6,1,2,'cm','jet');
if saveFigures
    set(gcf,'position',[0,0,400,400])
    saveas(gcf,'sec6-squ20-pt198','epsc')
end

%% c9: plot solution from secondary bpt3 close to film rupture
plotsol('square20-sec-bif3','pt238',3,1,2,'cm','jet');
if saveFigures
    set(gcf,'position',[0,0,400,400])
    saveas(gcf,'sec3-squ20-pt238','epsc')
end

%% c10: export of plot data
% extract and save solution data.

% load solution close to film-rupture
data_squ20_pt500 = loadp('squ20','pt500');
sol_squ20_pt500 = [data_squ20_pt500.pdeo.grid.p;data_squ20_pt500.u(1:end-3).'];

% load secondary bifurcation solution
data_squ20_6 = loadp('square20-sec-bif6','pt198');
sol_squ20_6_pt198 = [data_squ20_6.pdeo.grid.p;data_squ20_6.u(1:end-3).'];

data_squ20_3 = loadp('square20-sec-bif3','pt238');
sol_squ20_3_pt238 = [data_squ20_3.pdeo.grid.p;data_squ20_3.u(1:end-3).'];

% extract bifurcation branch data
data_init = loadpp('init-squ20');
brainit = data_init.branch([4,6],:);

data_squ20 = loadpp('squ20');
bra_squ20 = data_squ20.branch([4,6],:);

data_sec6_squ20 = loadpp('square20-sec-bif6');
bra_sec6_squ20 = data_sec6_squ20.branch([4,6],:);

data_sec3_squ20 = loadpp('square20-sec-bif3');
bra_sec3_squ20 = data_sec3_squ20.branch([4,6],:);

% branch data M = 8
data_init8 = loadpp('init-squ');
brainit8 = data_init8.branch([4,6],:);

data_squ = loadpp('squ');
bra_squ = data_squ.branch([4,6],:);

% write data to text files
if saveData
    formatSpec = '(%5.5f,%5.5f) ';

    fileID = fopen('branchDataSqu20.txt','w');
    fprintf(fileID,'%1s\n\n',['==================',' brainit ','==================']);
    fprintf(fileID,formatSpec,brainit);

    fprintf(fileID,'\n\n%1s\n\n',['==================',' bra_squ20 ','==================']);
    fprintf(fileID,formatSpec,bra_squ20);

    fprintf(fileID,'\n\n%1s\n\n',['==================',' sec6 ','==================']);
    fprintf(fileID,formatSpec,bra_sec6_squ20);

    fprintf(fileID,'\n\n%1s\n\n',['==================',' sec3 ','==================']);
    fprintf(fileID,formatSpec,bra_sec3_squ20);
    fclose(fileID);

    fileID = fopen('branchDataSquJoined.txt','w');
    fprintf(fileID,'%1s\n\n',['==================',' brainit20 ','==================']);
    fprintf(fileID,formatSpec,brainit);

    fprintf(fileID,'\n\n%1s\n\n',['==================',' bra_squ20 ','==================']);
    fprintf(fileID,formatSpec,bra_squ20);
    
    fprintf(fileID,'%1s\n\n',['==================',' brainit8 ','==================']);
    fprintf(fileID,formatSpec,brainit8);

    fprintf(fileID,'\n\n%1s\n\n',['==================',' bra_squ8 ','==================']);
    fprintf(fileID,formatSpec,bra_squ);
    fclose(fileID);

    fileID = fopen('solutionDataSqu20.txt','w');
    fprintf(fileID,'\n\n%1s\n\n',['==================',' solution squares20 pt500 ','==================']);
    fprintf(fileID,formatSpec,sol_squ20_pt500);
    fclose(fileID);

    
    fileID = fopen('solutionDataSqu20Sec.txt','w');
    fprintf(fileID,'%1s\n\n',['==================',' solution squares20 secondary bif 6 pt 198 ','==================']);
    fprintf(fileID,formatSpec,sol_squ20_6_pt198);

    fprintf(fileID,'\n\n%1s\n\n',['==================',' solution squares20 secondary bif 3 pt 238 ','==================']);
    fprintf(fileID,formatSpec,sol_squ20_3_pt238);

    fclose(fileID)
end