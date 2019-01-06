function [data] = initializer(const,v0,a0)
%% Tittle Block:

% Written by: Reed Clark
% Date Created: 1-4-2019
% Revised By: Reed Clark
% Revision Description: added v0 and a0 inputs for complete
% initialization

%% Description:
% This function initializes a structure for metadata storage

%% Inputs:
% const: structure, list of shared constants, controlled by the
% constants function

%% Outputs: 
% data: structure containing initialized metadata storage

%% Constants:
dt = const.dt;

fields = {'time' 'position' 'velocity' 'acceleration' 'alpha' 'beta' 'mass' 'hi_and_ho' 'hi_ratio_ho'};
subFields = {'fig' 'plotStyle' 'axes' 'xLabel' 'yLabel' 'value'};

%% Equations:
% Define equations used throughout code in LaTex for use in publishing.

%% Code:
t = 0:dt:30;

for i = 1:length(fields)
    if strcmpi(fields{i}, 'time') % must match fields definition!
        data.(fields{i}) = t;
    else
        for j = 1:length(subfields)
            switch subFields{j}
                case 'fig' % must match subFields definition!
                    data.(fields{i}).(subFields{j}) = figure('Name',fields{i},'NumberTitle','off');
                case 'plotStyle' % must match subFields definition!
                    data.(fields{i}).(subFields{j}) = '*';
                case 'axes' % must match subFields definition!
                    data.(fields{i}).(subFields{j}) = [0 inf 0 inf];
                case 'xLabel' % must match subFields definition!
                    data.(fields{i}).(subFields{j}) = 'time';
                case 'yLabel' % must match subFields definition!
                    data.(fields{i}).(subFields{j}) = fields{i};
                case 'value' % must match subFields definition!
                    data.(fields{i}).(subFields{j}) = zeros(1,length(t));
                    if srtcmpi(fields{i}, 'velocity') % must match fields definition!
                        data.(fields{i}).(subFields{j})(1) = v0;
                    elseif strcmpi(fields{i}, 'acceleration') % must match fields definition!
                        data.(fields{i}).(subFields{j})(1) = a0;
                    elseif strcmpi(fields{i}, 'hi_and_ho') % must match fields definition!
                        data.(fields{i}).(subFields{j}) = zeros(2,length(t));
                    end
            end
        end
    end
end
end