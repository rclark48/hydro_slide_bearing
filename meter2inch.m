function [out] = meter2inch(m)
%% Tittle Block:

% Written by: Reed Clark
% Date Created: 12-29-2018
% Revised By: Reed Clark
% Revision Description: initial commit
%% Description:
% Converts meters to inches

%% Inputs:
% in: can be a scalar, vector, or matrix

%% Outputs: 
% out: can be a scalar, vector, or matrix, dependent on the input

%% Constants:
conv = 25.4e-3; % [m --> in] conversion factor

%% Code:
out = m./conv;
end