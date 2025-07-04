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
p=gentau(p0,[1,1],'hex-up');        % Detected tangent directions are \phi_1 = cos(k_1*(x,y)) and \phi_2 = cos(k_2*(x,y))+cos(k_3*(x,y)). 
%                                     Generate tanget direction \phi_1 + \phi_2, which yields hexagons.
p.sol.ds=0.001;                     % positive for up-hexagons
p=pmcont(p,225);

%% c4: switch branch to periodic bifurcation branches and continue up to film-rupture through down-hexagons
p0=qswibra('init-hex','bpt1');
p=gentau(p0,[1,1],'hex-down');      
p.sol.ds=-0.001;                    % negative for down-hexagons
p=pmcont(p,276);

%% c5: plot bifurcation diagram hexagons
hold on;
plotbra('init-hex','cl','k');
plotbra('hex-up','cl','b');
plotbra('hex-down','cl','b');
hold off;
if saveFigures
    set(gcf,'position',[0,0,680,400])
    saveas(gcf,'bifurcation-diag-hex','epsc');
end

%% c6: plot up-hexagon close to the bifurcation point
plotsol('hex-up','pt5',1,1,2,'cm','jet');
if saveFigures
    set(gcf,'position',[0,0,680,400])
    saveas(gcf,'small-hex-up','epsc');
end

%% c7: plot up-hexagon at fold (M = 8.42325, pt50)
plotsol('hex-up','fpt1',1,1,2);
if saveFigures
    set(gcf,'position',[0,0,500,400])
    saveas(gcf,'fold-hex-up','epsc');
end

%% c8: plot up-hexagon close to film-rupture
plotsol('hex-up','pt224',1,1,3,'cm','jet');
if saveFigures
    set(gcf,'position',[0,0,680,400])
    saveas(gcf,'film-rupture-hex-up','epsc');
end

%% c9: plot down-hexagon close to the bifurcation point
plotsol('hex-down','pt5',1,1,2,'cm','jet');
if saveFigures
    set(gcf,'position',[0,0,680,400])
    saveas(gcf,'small-hex-down','epsc');
end

%% c10: plot down-hexagon close to film-rupture
plotsol('hex-down','pt275',1,1,3,'cm','jet');
if saveFigures
    set(gcf,'position',[0,0,680,400])
    saveas(gcf,'film-rupture-hex-down','epsc');
end

%% c11: plot integration constant \lambda = K(0) - K(v)
hold on;
plotbra('init-hex',15,3,'cl','k');
plotbra('hex-up',15,3,'cl','k');
plotbra('hex-down',15,3,'cl','k');
hold off;
if saveFigures
    set(gcf,'position',[0,0,680,400])
    saveas(gcf,'constant-K-hex','epsc');
end

%% c12: plot L^2 norm of log(1+v)
hold on;
plotbra('init-hex',16,5,'cl','k');
plotbra('hex-up',16,5,'cl','k');
plotbra('hex-down',16,5,'cl','k');
hold off;
if saveFigures
    set(gcf,'position',[0,0,680,400])
    saveas(gcf,'l2norm-hex','epsc');
end

%% c13: plot minimum
hold on;
plotbra('init-hex',17,4,'cl','k');
plotbra('hex-up',17,4,'cl','k');
plotbra('hex-down',17,4,'cl','k');
hold off;
if saveFigures
    set(gcf,'position',[0,0,680,400])
    saveas(gcf,'minimum-hex','epsc');
end

%% c14: export of plot data
% extract and save solution data.

% load solution up-hexagon close to the bifurcation point
data_hex_up_pt5 = loadp('hex-up','pt5');
sol_hex_up_pt5 = [data_hex_up_pt5.pdeo.grid.p;data_hex_up_pt5.u(1:end-3).'];

% load solution up-hexagon at fold
data_hex_up_fpt1 = loadp('hex-up','fpt1');
sol_hex_up_fpt1 = [data_hex_up_fpt1.pdeo.grid.p;data_hex_up_fpt1.u(1:end-3).'];

% load solution up-hexagon close to film-rupture
data_hex_up_pt224 = loadp('hex-up','pt224');
sol_hex_up_pt224 = [data_hex_up_pt224.pdeo.grid.p;data_hex_up_pt224.u(1:end-3).'];

% load solution down-hexagon close to the bifurcation point
data_hex_down_pt5 = loadp('hex-down','pt5');
sol_hex_down_pt5 = [data_hex_down_pt5.pdeo.grid.p;data_hex_down_pt5.u(1:end-3).'];

% load solution down-hexagon close to film-rupture
data_hex_down_pt275 = loadp('hex-down','pt275');
sol_hex_down_pt275 = [data_hex_down_pt275.pdeo.grid.p;data_hex_down_pt275.u(1:end-3).'];

% extract bifurcation branch data
data_init = loadpp('init-hex');
brainit = data_init.branch([4,6],:);

data_hex_up = loadpp('hex-up');
bra_hex_up = data_hex_up.branch([4,6],:);

data_hex_down = loadpp('hex-down');
bra_hex_down = data_hex_down.branch([4,6],:);

% extract other plots data
data_init_const_hex = loadpp('init-hex');
bra_init_const_hex = data_init_const_hex.branch([4,9],:);

data_init_l2norm_hex = loadpp('init-hex');
bra_init_l2norm_hex = data_init_l2norm_hex.branch([4,11],:);

data_init_min_hex = loadpp('init-hex');
bra_init_min_hex = data_init_min_hex.branch([4,10],:);


data_const_hex_up = loadpp('hex-up');
bra_const_hex_up = data_const_hex_up.branch([4,9],:);

data_l2norm_hex_up = loadpp('hex-up');
bra_l2norm_hex_up = data_l2norm_hex_up.branch([4,11],:);

data_min_hex_up = loadpp('hex-up');
bra_min_hex_up = data_min_hex_up.branch([4,10],:);


data_const_hex_down = loadpp('hex-down');
bra_const_hex_down = data_const_hex_down.branch([4,9],:);

data_l2norm_hex_down = loadpp('hex-down');
bra_l2norm_hex_down = data_l2norm_hex_down.branch([4,11],:);

data_min_hex_down = loadpp('hex-down');
bra_min_hex_down = data_min_hex_down.branch([4,10],:);

% write data to text files
if saveData
    formatSpec = '(%5.5f,%5.5f) ';

    fileID = fopen('branchDataHex.txt','w');
    fprintf(fileID,'%1s\n\n',['==================',' brainit ','==================']);
    fprintf(fileID,formatSpec,brainit);

    fprintf(fileID,'\n\n%1s\n\n',['==================',' bra_hex_up ','==================']);
    fprintf(fileID,formatSpec,bra_hex_up);

    fprintf(fileID,'\n\n%1s\n\n',['==================',' bra_hex_down ','==================']);
    fprintf(fileID,formatSpec,bra_hex_down);
    fclose(fileID);


    fileID = fopen('solutionDataHex.txt','w');
    fprintf(fileID,'%1s\n\n',['==================',' solution up-hexagon pt5 ','==================']);
    fprintf(fileID,formatSpec,sol_hex_up_pt5);

    fprintf(fileID,'\n\n%1s\n\n',['==================',' solution up-hexagon fp1 ','==================']);
    fprintf(fileID,formatSpec,sol_hex_up_fpt1);

    fprintf(fileID,'\n\n%1s\n\n',['==================',' solution up-hexagon pt224 ','==================']);
    fprintf(fileID,formatSpec,sol_hex_up_pt224);

    
    fprintf(fileID,'\n\n%1s\n\n',['==================',' solution down-hexagon pt5 ','==================']);
    fprintf(fileID,formatSpec,sol_hex_down_pt5);

    fprintf(fileID,'\n\n%1s\n\n',['==================',' solution down-hexagon pt275 ','==================']);
    fprintf(fileID,formatSpec,sol_hex_down_pt275);
    fclose(fileID);


    fileID = fopen('constantPlotHex.txt','w');
    fprintf(fileID,'%1s\n\n',['==================',' constant_init ','==================']);
    fprintf(fileID,formatSpec,bra_init_const_hex);

    fprintf(fileID,'\n\n%1s\n\n',['==================',' constant_up ','==================']);
    fprintf(fileID,formatSpec,bra_const_hex_up);

    fprintf(fileID,'\n\n%1s\n\n',['==================',' constant_down ','==================']);
    fprintf(fileID,formatSpec,bra_const_hex_down);
    
    
    fileID = fopen('l2normPlotHex.txt','w');
    fprintf(fileID,'%1s\n\n',['==================',' l2norm_init ','==================']);
    fprintf(fileID,formatSpec,bra_init_l2norm_hex);

    fprintf(fileID,'\n\n%1s\n\n',['==================',' l2norm_up ','==================']);
    fprintf(fileID,formatSpec,bra_l2norm_hex_up);

    fprintf(fileID,'\n\n%1s\n\n',['==================',' l2norm_down ','==================']);
    fprintf(fileID,formatSpec,bra_l2norm_hex_down);


    fileID = fopen('minHex.txt','w');
    fprintf(fileID,'%1s\n\n',['==================',' min_init ','==================']);
    fprintf(fileID,formatSpec,bra_init_min_hex);

    fprintf(fileID,'\n\n%1s\n\n',['==================',' min_up ','==================']);
    fprintf(fileID,formatSpec,bra_min_hex_up);

    fprintf(fileID,'\n\n%1s\n\n',['==================',' min_down ','==================']);
    fprintf(fileID,formatSpec,bra_min_hex_down);
end