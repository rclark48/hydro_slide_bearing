function [data] = initializer(const)
%% Tittle Block:

% Written by: Reed Clark
% Date Created: 1-4-2019
% Revised By: Reed Clark
% Revision Description: corrected year in date created

%% Description:
% This function initializes a structure for metadata storage

%% Inputs:
% const: structure, list of shared constants, controlled by the
% constants function

%% Outputs: 
% data: structure containing initialized metadata storage

%% Constants:
dt = const.dt;

fields = {'time' 'position' 'velocity' 'acceleration' 'alpha' 'beta'};
subFields = {'fig' 'plotStyle' 'axes' 'xLabel' 'yLabel' 'valVec'};

%% Equations:
% Define equations used throughout code in LaTex for use in publishing.

%% Code:
t = 0:dt:30;

for i = 1:length(fields)
    if i == 1
        data.(fields{i}) = t;
    else
        for j = 1:length(subfields)
            switch j
                case 1
                    data.(fields{i}).(subFields{j}) = figure('Name',fields{i},'NumberTitle','off');
                case 2
                    data.(fields{i}).(subFields{j}) = '*';
                case 3
                    data.(fields{i}).(subFields{j}) = [0 inf 0 inf];
                case 4
                    data.(fields{i}).(subFields{j}) = 'time';
                case 5
                    data.(fields{i}).(subFields{j}) = fields{i};
                case 6
                    data.(fields{i}).(subFields{j}) = zeros(1,length(t));
            end
        end
    end
end
end