function [data] = downTheRamp(filename)
%% Tittle Block:

% Written by: Reed Clark
% Date Created: 1-3-2019
% Revised By: Reed Clark
% Revision Description: added all relevant function calls and started
% building test algorithm.

%% Description:
% Given initial conditions, this function solves force and moment
% equilibrium equations as an unconstrained hydrodynamic slider
% bearing moves down an inclined pane. At each time step, it will
% solve for the inlet and outlet lubricant heights and the values or
% alpha and beta (representing the volume of material to be removed to
% properly locate the bearing center of mass). Once the bearing has
% traveled all the way down the ramp, the average value of alpha and
% beta will be chosen and then the bearing motion will be recalculated
% with a fixed value of alpha and beta (since the dimensions are
% technically invariant once machined in to the block). This method is
% the first attempt at finding the optimal bearing geometry.

%% Inputs:
% filename: string, specify a file name to save data to for post
% processing

%% Outputs: 
% 

%% Constants:
const = constants();

[a0, v0] = initialConditions(const);

disc = discretizer(25,const);

data = initializer(const,v0,a0);
%% Equations:
% Define equations used throughout code in LaTex for use in publishing.

%% Code:
step = 1;
iter = 1;

% Terminating condition initialization:
s = data.position.value(1);
a = 1;

% H initialization:
ho = 50; % [microns]
hi = 150; % [microns]

% Removed material initialization:
alpha = .5; % [inches]
beta = .25; % [inches]

while s >= -inch2meter(24) && a >= const.err.resid
    [pressure,H] = pressureSolver(U,H,disc,const,solver);
    com = COM(alpha, beta, const);
    [force, moment] = integrator(pressure,U,disc,H,const,com);
    res = residuals(force,moment,com,const);
    
    [data, res] = kinematics(step,data,const,force,com,res);
end
end