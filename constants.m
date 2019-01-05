function [const] = constants()
%% Tittle Block:

% Written by: Reed Clark
% Date Created: 12-30-2018
% Revised By: Reed Clark
% Revision Description: added residual tolerance and time step

%% Description:
% This is a modularization function that initializes all constants
% that are shared between multiple functions--ideally things such as
% design problem constraints.

%% Inputs:
% There are no explicit inputs to this function

%% Outputs: 
% const: structure, list of shared constants, controlled by the
% constants function

%% Constants:
const.A = inch2meter(1); % [in --> m]
const.B = inch2meter(2); % [in --> m]
const.L = inch2meter(2); % [in --> m]

const.g = 9.81; % [m/s^2]
const.theta = 2.6; % [degrees]
const.z0 = 2e-3; % [m]

const.density = 7750; %[kg/m^3], mild steel
const.viscosity = 54.3e-3; %[Pa-s], soybean oil ~75 F

% Relative percent error tolerance for Gauss-Seidel pressure solver
const.err.gauss = eps;
% Residual error tolerance for force and moment balance equations. The
% residuals are expressed as a normalized approximate percent error
% percentage, so the tolerance corresponds to roughly 1/10% error.
const.err.resid = 1e-3; 

% For most human timed events (such as this trial) are probably only
% accurate to a hundreth of a second, so running the simulation at one
% higher degree of resolution in the time domain should ensure that a
% estimated travel time truncated to the hundreth of a second will be
% a reasonable prediction.
const.dt = 1e-3; % [s]

%% Equations:
% There are no formal equations to define in this function.

%% Code:
% There is no explicit code in this function.
end