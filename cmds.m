%% main command file for the bifurcation analysis of a ---2D--- thermocapillary thin-film equation
% 
% The command file performs a numerical continuation of the second order
% equation
%   0 = \Delta h - (g(1+h)-(3/2)*M*(log((1+h)/(2+h))+1/(2+h))+K(0)+lambda
% Here, h is the fluid height, g is the gravitation constant, M is the
% Marangoni number (the primary bifurcation parameter), lambda is a
% constaint to conserve the mass of the bifurcating solution.
%%%%%%%%%%%%%%%%%%%% NOTE K and h here are -MK and v in the analysis 
%%%%%%%%%%%%%%%%%%%% TYPO h should be v
%
% The numerical continuation is performed with pde2path with Neuman
% boundary conditions on a domain (-lx,lx)x(-ly,ly).


%% clear mainspace (keep paths for pde2path) and close all plot windows
clc;
close all;
keep pphome;

% set to 1 to save plotted figures as eps
saveFigures = 0;
% set to 1 to save plot data
saveData = 0;

%% c1: init and set parameters
p = [];
lx = [7*pi 7*pi];   % set domain size
nx = [30 30];     % number of discretisation points per dimension in domain
Minit = 0;          % initial Marangoni number
ginit = 1;          % set gravitiy constant
lambdaInit = 0;     % set initial integration constant (lambda = K - K(0))
par = [Minit, ginit, lambdaInit];
p = tfinit(p,lx,nx,par);
p = setfn(p,'init');
para = 1;           % set Marangoni number as bifurcation parameter
p.nc.dsmax=0.03;%%%%%%%%%%%%%%%%%%%% copied
                                                                                                    
%% c2: continuation of the trivial branch
p = cont(p,200);%%%%%%%%%%%%%%%%%%%% 90 not enough

%% c3: switch branch to periodic bifurcation branches and continue 
p = swibra('init','bpt2','1D2'); p = cont(p,20); % branch at M^*(k_0)

p = swibra('init','bpt4','1D4'); p = cont(p,20); % branch at M^*(2k_0)

p = swibra('init','bpt6','1D6'); p = cont(p,20); % branch at M^*(3k_0)


%% c6: plot bifurcation diagram
hold on;
plotbra('init','cl','k');
plotbra('1D2','cl','b','lab',[1,10,15,16,17,18,19,20]);
plotbra('1D4','cl','b','lab',[5,15]);%%%%%%%%%%%%%%%%%%%% was ms instead of lab
plotbra('1D6','cl','b','lab',[5,20]);%%%%%%%%%%%%%%%%%%%% was ms instead of lab
hold off;

% % save bifurcation diagram as eps
% if saveFigures
%     set(gcf,'position',[0,0,500,400])
%     saveas(gcf,'bifurcation-graph','epsc');
% end

%% c7 : plot solutions
% % plot bifurcating solution close to the bifurcation point
p2 = loadp('1D2','pt5');
% figure(4) %%%%%%%%%%%%%%%%%%%% commented out
% plot solution shifted by a half-period (so that minimum is at x=0)
templx=lx(1)
plot(p2.pdeo.grid.p+templx*ones(size(p2.pdeo.grid.p)),p2.u(1:end-3),'k',...
    p2.pdeo.grid.p-templx*ones(size(p2.pdeo.grid.p)),p2.u(1:end-3),'k');
xlim([-lx(1),lx(1)]);
ylim([-lx(2),lx(2)]);
title('Solution plot at pt5');

% % save plot as eps
% if saveFigures
%     set(gcf,'position',[0,0,500,200])
%     saveas(gcf,'sol-pt5','epsc');
% end


% % plot 'film-rupture' solution
% p2 = loadp('1D2','pt40');
% figure(5)
% % plot solution shifted by a half-period (so that minimum is at x=0)
% plot(p2.pdeo.grid.p+lx*ones(size(p2.pdeo.grid.p)),p2.u(1:end-3),'k',...
%     p2.pdeo.grid.p-lx*ones(size(p2.pdeo.grid.p)),p2.u(1:end-3),'k');
% xlim([-lx,lx]);
% ylim([-1,0.5]);
% title('Solution plot at pt40');

% % save plot as eps
% if saveFigures
%     set(gcf,'position',[0,0,500,200])
%     saveas(gcf,'sol-pt40','epsc');
% end

% %% export of plot data
% % extract and save solution data.
% 
% % interpolation grid with finer mesh close to 0 (where the minimum will be located in the plot).
% % This allows to reduce the number of points necessary for the plot to look 'smooth'
% xInterpLeft = [linspace(-lx,-0.5,50),linspace(-0.5,0,200)];
% xInterpRight = [linspace(0,0.5,200),linspace(0.5,lx,50)];
% 
% % set up temporary fine mesh
% xTempLeft = linspace(-lx,0,1000);
% xTempRight = linspace(0,lx,1000);
% 
% % load solution close to bifurcation point
% data_1D2_pt5 = loadp('1D2','pt5');
% % interpolate solution to temporary mesh and shifting solution by a half-period
% soltemp = [interp1(data_1D2_pt5.pdeo.grid.p,data_1D2_pt5.u(1:end-3),xTempLeft(1:end-1)),interp1(data_1D2_pt5.pdeo.grid.p,data_1D2_pt5.u(1:end-3),xTempRight(2:end))];
% % interpolate to final mesh and save grid points and solution points
% solpt5 = [xInterpLeft,xInterpRight;...
%     interp1([xTempRight(2:end),xTempLeft(1:end-1)],soltemp,xInterpLeft),interp1([xTempRight(2:end),xTempLeft(1:end-1)],soltemp,xInterpRight)];
% 
% % load 'film-rupture' solution
% data_1D2_pt40 = loadp('1D2','pt40');
% % interpolate solution to temporary mesh and shifting solution by a
% % half-period
% soltemp = [interp1(data_1D2_pt40.pdeo.grid.p,data_1D2_pt40.u(1:end-3),xTempLeft(1:end-1)),interp1(data_1D2_pt40.pdeo.grid.p,data_1D2_pt40.u(1:end-3),xTempRight(2:end))];
% % interpolate to final mesh and save grid points and solution points
% solpt40 = [xInterpLeft,xInterpRight;...
%     interp1([xTempRight(2:end),xTempLeft(1:end-1)],soltemp,xInterpLeft),interp1([xTempRight(2:end),xTempLeft(1:end-1)],soltemp,xInterpRight)];
% 
% 
% % extract bifurcation branch data
% data_init = loadpp('init');
% brainit = data_init.branch([4,6],:);
% 
% data_1D2 = loadpp('1D2');
% bra1D2 = data_1D2.branch([4,6],:);
% 
% data_1D4 = loadpp('1D4');
% bra1D4 = data_1D4.branch([4,6],:);
% 
% data_1D6 = loadpp('1D6');
% bra1D6 = data_1D6.branch([4,6],:);
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
%     fprintf(fileID,formatSpec,bra1D2);
% 
%     fprintf(fileID,'\n\n%1s\n\n',['==================',' bra1D4 ','==================']);
%     fprintf(fileID,formatSpec,bra1D4);
% 
%     fprintf(fileID,'\n\n%1s\n\n',['==================',' bra1D6 ','==================']);
%     fprintf(fileID,formatSpec,bra1D6);
%     fclose(fileID);
% 
%     fileID = fopen('solutionData.txt','w');
%     fprintf(fileID,'%1s\n\n',['==================',' solution pt5 ','==================']);
%     fprintf(fileID,formatSpec,solpt5);
% 
%     fprintf(fileID,'\n\n%1s\n\n',['==================',' solution pt40 ','==================']);
%     fprintf(fileID,formatSpec,solpt40);
%     fclose(fileID);
% end