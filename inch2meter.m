function [m] = inch2meter(in)
%% Title Block (update when revised)
%            Written by: Reed Clark
%          Date Created: 12-18-2018
%            Revised By: Reed Clark
%  Revision Description: title block updates: removed lines
%  about revision date and revision number--git tracks this,
%  so no need to manually update
%% Description:
% Convert units from inches to meters.

%% Inputs:
% in: can be a scalar, vector, or matrix

%% Code:
m = in*25.4e-3;
end