function [const] = constants()
%% Tittle Block:

% Written by: Reed Clark
% Date Created: 12-30-2018
% Revised By: Reed Clark
% Revision Description: corrected output description

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

const.err.gauss = eps; % relative error tolerance for Gauss-Seidel pressure solver

%% Equations:
% There are no formal equations to define in this function.

%% Code:
% There is no explicit code in this function.
end