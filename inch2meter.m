function [out] = inch2meter(in)
%% Tittle Block:

% Written by: Reed Clark
% Date Created: 12-29-2018
% Revised By: Reed Clark
% Revision Description: "save-as" to get formatting up to snuff with
% template

%% Description:
% Converts inches to meters

%% Inputs:
% in: can be a scalar, vector, or matrix

%% Outputs: 
% out: can be a scalar, vector, or matrix, dependent on the input

%% Constants:
conv = 25.4e-3; % [in --> m] conversion factor

%% Code:
out = in*conv;
end