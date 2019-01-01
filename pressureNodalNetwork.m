function [pressure,coord] = pressureNodalNetwork(coord,height,const,solver)
%% Tittle Block:

% Written by: Reed Clark
% Date Created: 12-31-2018
% Revised By: Reed Clark
% Revision Description: initial commit

%% Description:
% The hydrodynamic pressure is given by the discretized Reynold's
% Equation in 2 dimensions. This function defines the discretization
% and solves for the pressure distribution. The function can solve for
% the nodal pressures in two different ways: one using Gauss-Seidel
% iteration and another just using mldivide and letting MATLAB decide
% how to solve the pressure distribution. Regardless of the solution
% method selected, the function will return the result formatted in
% the same way. The idea of providing two options is to run some test
% cases and see which version of the solver runs faster. The system
% should be diagonally dominant and therefore have a guaranteed
% convergent solution, but the Gauss-Seidel solver may be neccessary
% as the resolution (i.e. the number of nodal pressure equations)
% increases. There is also the third option of using both solvers and
% then comparing the two resluts as a sanity check.

%% Inputs:
% coord: structure, fields define the resolution of the nodal network 

% height: structure, fields define the inlet height, outlet height,
% slope of bearing converging gap, and the ratio of inlet to outlet
% heights

% const: structure, list of shared constants, controlled by the
% constants function

% solver: string, defines which solver to use: Gauss-Seidel, mldivide,
% or both as a sanity check

%% Outputs: 
% pressure: matrix <m x n>, [Pa], nodal pressure distribution

% coord.discretized: structure.field, cell array <m x n>, [m], each
% bucket contains a vector <1 x 3> with the x, y, and z (hydrodynamic
% height) nodal coordinates as disctretized by the specified
% resolutiuon

%% Constants:
B = const.B;
L = const.L;

%% Equations:
% Define equations used throughout code in LaTex for use in publishing.

% Block equation:
%% 
% 
% $$e^{\pi i} + 1 = 0$$
% 

% Inline equation:
%%
% $x^2+e^{\pi i}$

%% Code:

end