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
Minit = 4.9;                        % initial Marangoni number
ginit = 1;                          % set gravitational constant
lambdaInit = 0;                     % set initial integration constant (lambda = M*( K(0) - K ))
par = [Minit, ginit, lambdaInit];
p = tfinit(p,lx,nx,par);
p = setfn(p,'init-squ');
para = 1;                           % set Marangoni number as bifurcation parameter
p.nc.dsmax=0.03;
                                                                                                    
%% c2: continuation of the trivial branch
p = cont(p,15);                     % set up to detect bifurcation point at M = 5 (one full square)

%% c3: switch branch to periodic bifurcation branches and continue up to film-rupture through squares
p0=cswibra('init-squ','bpt1'); 
p=seltau(p0, 2,'squ',3); 
p.sol.ds=0.01; 
p=pmcont(p,120); 

%% c-sec: secondary bifurcation
p = swibra('squ','bpt2','square-sec-bif');
p.sol.ds=0.01;
p=pmcont(p,100);

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

%% c7: plot secondary bifurcations
if saveFigures
    set(gcf,'position',[0,0,500,400])
    plotsol('square-sec-bif','pt100',1,1,2,'cm','jet')
end
saveas(gcf,'squ-secondary-bif','epsc')

%% c7: export of plot data
% extract and save solution data.

% interpolation grid with finer mesh close to ---corners--- (where the minimum will be located in the plot).
% This allows to reduce the number of points necessary for the plot to look 'smooth'
% xmax=lx(1);
% ymax=lx(2);

% xInterpLeft = [linspace(-xmax,-xmax+0.5,200),linspace(-xmax+0.5,0,50)];
% xInterpRight = [linspace(0,xmax-0.5,50),linspace(xmax-0.5,xmax,200)];
% yInterpBott = [linspace(-ymax,-ymax+0.5,200),linspace(-ymax+0.5,0,50)];
% yInterpTop = [linspace(0,ymax-0.5,50),linspace(ymax-0.5,ymax,200)];
% 
% % set up temporary fine mesh
% xTempLeft = linspace(-xmax,0,1000);
% xTempRight = linspace(0,xmax,1000);
% yTempBott = linspace(-ymax,0,1000);
% yTempTop = linspace(0,ymax,1000);

% load solution close to bifurcation point
% data_squ_pt5 = loadp('squ','pt5');
% interpolate solution to temporary mesh
% soltemp = [interp2(data_squ_pt5.pdeo.grid.p(1,:),data_squ_pt5.pdeo.grid.p(2,:),data_squ_pt5.u(1:end-3),xTempLeft(1:end-1),yTempBott(1:end-1)),...   % bottom left
%            interp2(data_squ_pt5.pdeo.grid.p(1,:),data_squ_pt5.pdeo.grid.p(2,:),data_squ_pt5.u(1:end-3),xTempRight(2:end),yTempBott(1:end-1)),...    % bottom right
%            interp2(data_squ_pt5.pdeo.grid.p(1,:),data_squ_pt5.pdeo.grid.p(2,:),data_squ_pt5.u(1:end-3),xTempLeft(1:end-1),yTempTop(2:end)),...      % top left
%            interp2(data_squ_pt5.pdeo.grid.p(1,:),data_squ_pt5.pdeo.grid.p(2,:),data_squ_pt5.u(1:end-3),xTempRight(2:end),yTempUp(2:end))];          % top right
% interpolate to final mesh and save grid points and solution points
% sol_squ_pt5 = [xInterpLeft,xInterpRight,yInterpBott,yInterpTop;...
%     interp2([xTempLeft(1:end-1),xTempRight(2:end),yTempBott(1:end-1),yTempTop(2:end)],soltemp,xInterpLeft,yInterpBott),...
%     interp2([xTempLeft(1:end-1),xTempRight(2:end),yTempBott(1:end-1),yTempTop(2:end)],soltemp,xInterpRight,yInterpBott),...
%     interp2([xTempLeft(1:end-1),xTempRight(2:end),yTempBott(1:end-1),yTempTop(2:end)],soltemp,xInterpLeft,yInterpTop),...
%     interp2([xTempLeft(1:end-1),xTempRight(2:end),yTempBott(1:end-1),yTempTop(2:end)],soltemp,xInterpRight,yInterpTop)];

% [X,Y] = meshgrid(linspace(-2*pi,2*pi,500),linspace(-2*pi,2*pi,500));
% Z = interp2(data_squ_pt5.pdeo.grid.p(1,:),data_squ_pt5.pdeo.grid.p(2,:),data_squ_pt5.u(1:end-3),X,Y);
% 
% fig = figure('units','pixels','position',[100 200 1000 1000]);
% axes(fig,'units','pixels','position',[0 0 1000 1000]);
% 
% [X,Y] = meshgrid(data_squ_pt5.pdeo.grid.p(1,:),data_squ_pt5.pdeo.grid.p(2,:));
% 
% s = surfc(data_squ_pt5.pdeo.grid.p(1,:),data_squ_pt5.pdeo.grid.p(2,:),data_squ_pt5.u(1:end-3));
% axis off;
% view(2);
% pbaspect([xmax/ymax 1 1])
% grid off;
% s(1).EdgeColor="none";
% colormap jet
% % clim([-maxAmp,maxAmp]);
% ax = gca;
% ax.ZLim(2) = max(f(X,Y),[],"all");
% s(2).ZLocation = 'zmax';
% s(2).LevelList = heights;
% s(2).EdgeColor="k";


% % extract bifurcation branch data
% data_init = loadpp('init-squ');
% brainit = data_init.branch([4,6],:);
% 
% data_squ = loadpp('squ');
% bra_squ = data_1D2.branch([4,6],:);
% 
% % write data to text files
% if saveData
%     formatSpec = '(%5.5f,%5.5f) ';
% 
%     fileID = fopen('branchData.txt','w');
%     fprintf(fileID,'%1s\n\n',['==================',' braInit ','==================']);
%     fprintf(fileID,formatSpec,brainit);
% 
%     fprintf(fileID,'\n\n%1s\n\n',['==================',' bra1D2 ','==================']);
%     fprintf(fileID,formatSpec,bra_squ);
%     fclose(fileID);
% 
%     fileID = fopen('solutionData.txt','w');
%     fprintf(fileID,'%1s\n\n',['==================',' solution squares pt5 ','==================']);
%     fprintf(fileID,formatSpec,sol_squ_pt5);
%     fclose(fileID);
% end