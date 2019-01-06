function [data, res] = kinematics(step,data,const,force,com,res)
%% Tittle Block:

% Written by: Reed Clark
% Date Created: 1-5-2019
% Revised By: Reed Clark
% Revision Description: initial commit

%% Description:
% This function computes the kinmatics of the bearing moving down the
% ramp at each time step

%% Inputs:
% step: scalar, [unitless],index of the current time step

% data: structure containing initialized metadata storage

% const: structure, list of shared constants, controlled by the
% constants function

% force.drag: scalar, [N], hydrodynamic drag

% com.mass: scalar, [kg], mass of bearing

%% Outputs: 
% data: structure containing initialized metadata storage

%% Constants:
a0 = data.acceleraction.value(1);
v = data.velocity.value(step);
s = data.position.value(step);

dt = const.dt;

D = force.drag;

m = com.mass;

%% Equations:
% Acceleration:
% a_{t+1} = \frac{D_t}{m_t} + a_0

% Velocity:
% v_{t+1} = v_t + a_{t+1} dt

% Position:
% s_{t+1} = s_t + v_t dt + \frac{1}{2} a_{t+1} dt^2

%% Code:
% Acceleration:
% a0 is negative and will yield aPlus as negative, which is consistent
% with the assigned coordinate system
aPlus = D./m + a0;
data.acceleration.value(step + 1) = aPlus;
res.acceleration = abs(aPlus./a0);

% Velocity:
% Based on the assigned coordinate system, this will be an
% increasingly negative value that should eventually level off as
% steady state is approached.
vPlus = v + aplus*dt;
data.velocity.value(step + 1) = vPlus;

% Position:
% Based on the assigned coordinate system, this will be an
% increasingly negative value--eventually reaching -24 inches when it
% reaches the end of the ramp.
sPlus = s + v*dt + 1/2*aPlus*dt^2;
data.position.value(step + 1) = sPlus;

end