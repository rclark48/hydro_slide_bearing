function [m] = inch2meter(in)
%% Title Block (update when revised)
%            Written by: Reed Clark
%          Date Created: 12-18-2018
%       Revision Number: 0
% Date of Last Revision: 12-24-2018
%            Revised By: Reed Clark
%  Revision Description: modified title block for git source
%  control

%% Description:
% Convert units from inches to meters.

%% Inputs:
% in: can be a scalar, vector, or matrix

%% Code:
m = in*25.4e-3;
end