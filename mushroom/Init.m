%% Initialize a flat membrane patch for use with the axisymmetric Helfrich membrane model


% Inputs:
%   alpha - dimensionless patch area
%   mesh - meshing for the domain, runs from 0 to 1, i.e. 0:0.01:1
%   lambda - membrane tension at the boundary, in units of pN/nm
%   k0 - bending rigidity of bare membrane, in units of pN*nm
%   R0 - nondimensionalization length

% Output:
%   initSol - initialized solution array


function initSol = Init(alpha, lambda, k0, R0)
mesh=(0:0.00005:1).^4;

ds=0.0001;    % Small deviation from x=0 (necessary to avoid division by 0)

t = alpha.*mesh;     % area mesh points

lam = lambda*R0^2/k0;   % dimensionless membrane tension

initSol = [ds + sqrt(2*t)           % x
           zeros(1, length(t))      % y
           zeros(1, length(t))      % psi
           zeros(1, length(t))      % h
           zeros(1, length(t))      % l  
           lam*ones(1, length(t))]; % lambda