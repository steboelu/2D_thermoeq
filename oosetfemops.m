function p = oosetfemops(p)
% OOSETFEMOPS set up function to save FEM mass and stiffness matrix in
% state variable p
%
%   p = oosetfemops(p)
%
% Input:
%   p : 'p' variable of pde2path (contains collection of current
%       continuation data)
%
% Ouput:
%   p : updated 'p' variable of pde2path
%

gr=p.pdeo.grid;%%%%%%%%%%%%%%%%%%%% copy data before using
[K,M,~] = p.pdeo.fem.assema(gr,1,1,1); % assemble matrices
p.mat.vM=sum(M,1); % such that vM*u=\int u dx
p.mat.M = M; % mass matrix
p.mat.K = K; % stiffness matrix
%%%%%%%%%%%%%%%%%%%% unclear to me what mass matrix is

%%%%%%%%%%%%%%%%%%%% copied but probably not needed
% [Dx,Dy]=p.pdeo.fem.gradientMatrices(gr); 
% c2p=center2PointMatrix(gr); 
% 
% Dx=c2p*Dx; 
% Dy=c2p*Dy; 
% Dz=0*Dy;
% 
% p.mat.Dx=Dx; 
% p.mat.Dy=Dy; 
% p.mat.Dz=Dz;
%%%%%%%%%%%%%%%%%%%%
end