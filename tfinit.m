function p = tfinit(p,lx,nx,par)
% TFINIT initializes pde2path parameters in thin-film (tf) problem
%
%   p = tfinit(p,lx,nx,par)
%
% Input:
%   p : 'p' variable of pde2path (contains collection of current
%       continuation data)
%   lx: to define the domain (-lx(1),lx(1)) x (-lx(2),lx(2))
%   nx: number of (spatial) grid points per dimension             
%   par : initial values of problem parameters (here: M,g,lambda)
%
% Output:
%   p : initialized 'p' variable of pde2path
%

p = stanparam(p);                       % assign 'standard values' to p
screenlayout(p);                        % set standard positions for pde2path figures
p.nc.neq = 1;                           % set number of equations
p.sw.sfem = -1;                         % flags use of oopde
p.plot.auxdict = {'M','g','lambda','min(v)','L^2 log(1+v)'};    % parameter names and additional plots

p.fuha.sG = @sG;                        % set bifurcation function
p.fuha.sGjac = @sGjac;                  % set derivative of bifurcation function
p.sw.jac=1;                             % 1 for analytic jacobian

sw.sym=2;                                       % criss-cross mesh style
pde = stanpdeo2D(lx(1),lx(2),nx(1),nx(2),sw);   % create 2D pde-object with (1D) mesh size h = 2*lx/nx

p.plot.pstyle=2;                        % contour plot
p.dim=2; 

p.vol = 4*lx(1)*lx(2);                  % set domain volume
p.pdeo = pde;                           % write pde-object into main struct object p
p.np = pde.grid.nPoints;                % number of grid points
p.nu = p.np*p.nc.neq;                   % number of nodal values
p.sol.xi = 1/p.nu;                      % assign norm weight (for internal pde2path use)

u = zeros(p.np,1);                      % initial (trivial) solution
p.u = [u;par'];                         % join solution and parameters in one vector

p = oosetfemops(p);                     % generate FEM matrices and assign to p
p.nc.ilam = [1,3];                      % set active parameters: M is primary bifurcation parameter and lambda is contraint parameter to obtain mass conservation
p.sol.ds = 0.01;                        % set initial arclength step size for continuation
p.nc.dsmax = 0.1;                       % set maximal arclength step size
p.nc.lammin = 0;                        % minimal bound for primary bifurcation parameter M
p.nc.neig = 30;                         % number of eigenvalues to compute
p.sw.foldcheck = 1;                     % turn on detection of foldpoints
p.sw.bifcheck = 2;                      % turn on bifurcation detection via counting eigenvalues
p.nc.nsteps = 100;                      % set default number of continuation steps

% add conservation of mass
p.nc.nq = 1;                            % number of constraints
p.fuha.qf = @qf;                        % set mass function
p.fuha.qfder = @qfder;                  % set derivative of mass function
p.sw.qjac=1;                            % 1 for analytic jacobian
end