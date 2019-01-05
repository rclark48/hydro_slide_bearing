function [res] = residuals(force,moment,com,const)
%% Tittle Block:

% Written by: Reed Clark
% Date Created: 1-5-2019
% Revised By: Reed Clark
% Revision Description: initial commit

%% Description:
% Template for uniform functions.

%% Inputs:
% force.loadSupport: scalar, [N], value of hydrodynamic load support

% force.drag: scalar, [N], hydrodynamic drag

% moment.loadSupport: scalar, [Nm], moment about com due to
% hydrodynamic load support.

% com.z: scalar, [m], z-coordinate of the bearing center of mass

% com.mass: scalar, [kg], mass of bearing

% const: structure, list of shared constants, controlled by the
% constants function

%% Outputs: 
% res: structure, [unitless], normalized residuals

%% Constants:
F = force.loadSupport;
D = force.drag;
M = moment.loadSupport;

z_cm = com.z;
m = com.mass;

g = const.g;
theta = const.theta;

%% Equations:
% Define equations used throughout code in LaTex for use in publishing.

% Sum of forces in zr direction:
% 0 = F - m g \cos \theta

% Sum of moments about COM:
% 0 = \int_{0}^{L} \int_{0}^{x_{cm}} p(x,y) x dx dy - \int_{0}^{L} \int_{x_{cm}}^{B} p(x,y) x dx dy - D z_{cm}

%% Code:
res.forceZ = (F - m*g*cosd(theta))./(m*g*cosd(theta));
res.momentCM = (M - D*z_cm)./M;
end