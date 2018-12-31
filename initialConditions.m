function [a0, v0] = initialConditions()
%% Tittle Block:

% Written by: Reed Clark
% Date Created: 12-30-2018
% Revised By: Reed Clark
% Revision Description: initial commit

%% Description:
% Computes the initial velocity of the bearing from when it's dropped
% to when it hits the ramp

%% Inputs:
% There are no inputs to this function--it is defined in terms of
% constants, defined further down in the constants block

%% Outputs: 
% a0: scalar, [m/s^2], initial acceleration of the bearing in the -x_r
% direction (the direction moving down the ramp in the ramp coordinate
% system)

% v0: scalar, [m/s], initial velocity of the bearing in the -x_r
% direction (the direction moving down the ramop in the ramp
% coordinate system)

%% Constants:
g = 9.81; % [m/s^2], gravitational acceleration
theta = 2.6; % [degrees], ramp inclination wrt to inertial frame
z0 = 2e-3; % [m], initial drop height

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