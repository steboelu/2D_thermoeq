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
nx = [400 400];                     % number of discretisation points per dimension in domain
Minit = 19.9;                       % initial Marangoni number
ginit = 1;                          % set gravitational constant
lambdaInit = 0;                     % set initial integration constant (lambda = M*( K(0) - K ))
par = [Minit, ginit, lambdaInit];   
p = tfinit(p,lx,nx,par);
p.fuha.outfu = @tfbra;
p = setfn(p,'init-hex20down');
para = 1;                           % set Marangoni number as bifurcation parameter
p.nc.dsmax=0.03;
                                                                                                    
%% c2: continuation of the trivial branch
p = cont(p,15);                    

%% c3: switch branch to periodic bifurcation branches and continue up to film-rupture through down-hexagons
p0=qswibra('init-hex20down','bpt1');
p=gentau(p0,[-1,1],'hex20down'); 
p.sol.ds=-0.001;                     
p=pmcont(p,400);        

% secondary bifurcations: no admissible patterns found