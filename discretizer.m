function [disc] = discretizer(n,const)
%% Tittle Block:

% Written by: Reed Clark
% Date Created: 1-3-2018
% Revised By: Reed Clark
% Revision Description: initial commit

%% Description:
% This function defines the grid discretization.

%% Inputs:
% n: scalar, fields define the resolution of the nodal network, i.e.
% the step size of the discretization

% const: structure, list of shared constants, controlled by the
% constants function


%% Outputs: 
% disc: structure, [m], nodal coordinates as disctretized by the
% specified resolutiuon. for code compatness, vectors and matrices
% will be used during computations while the coord array is just used
% to compatly pass information relevant to other functions out of the
% function

%% Constants:
B = const.B;
L = const.L;

%% Equations:
% Define equations used throughout code in LaTex for use in publishing.

%% Code:
% Grid discretization:
x = linspace(0,B,n); % vector of x-coordinates
dx = x(2)-x(1);
    disc.n = n;
    disc.x = x;
    disc.dx = dx;

m = round(L/B*n); % number of y-direction discretizations
y = linspace(0,L,m); % vector of y-coordinates
dy = y(2)-y(1);
    disc.m = m;
    disc.y = y;
    disc.dy = dy;

lambda = dx/dy;
    disc.lambda = lambda;

[x_mn, y_mn] = meshgrid(x,y);
    disc.x_mn = x_mn;
    disc.y_mn = y_mn;
end