function [a0, v0] = initialConditions(const)
%% Tittle Block:

% Written by: Reed Clark
% Date Created: 12-30-2018
% Revised By: Reed Clark
% Revision Description: added const input and redefined the constants

%% Description:
% Computes the initial velocity of the bearing from when it's dropped
% to when it hits the ramp

%% Inputs:
% const: structure, list of shared constants, controlled by the
% constants function

%% Outputs: 
% a0: scalar, [m/s^2], initial acceleration of the bearing in the -x_r
% direction (the direction moving down the ramp in the ramp coordinate
% system)

% v0: scalar, [m/s], initial velocity of the bearing in the -x_r
% direction (the direction moving down the ramp in the ramp
% coordinate system)

%% Constants:
g = const.g;
theta = const.theta;
z0 = const.z0;

%% Equations:
% Define equations used throughout code in LaTex for use in publishing.

% a0:
% a_0 = g\sin{\theta}

% v0:
% v_0 = g\sin{\theta}\sqrt{\frac{2z_0}{g\cos{\theta}}}

%% Code:
% a0:
a0 = g*sind(theta);

% v0:
v0 = g*sind(theta)*sqrt(2*z0/(g*cosd(theta)));
end