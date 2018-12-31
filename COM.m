function [x_cm,z_cm,m] = COM(alpha, beta, const)
%% Tittle Block:

% Written by: Reed Clark
% Date Created: 12-30-2018
% Revised By: Reed Clark
% Revision Description: added const input and redefined the constants

%% Description:
% Computes the bearing center of mass coordinates 

%% Inputs:
% alpha: scalar, [in]

% beta: scalar, [in]

% const: structure, list of shared constants, controlled by the
% constants function

alpha = inch2meter(alpha); % [in --> m]
beta = inch2meter(beta); % [in --> m]

%% Outputs: 
% x_cm: scalar, [m], x-coordinate of the bearing center of mass

% z_cm: scalar, [m], z-coordinate of the bearing center of mass

% mass: scalar, [kg], mass of bearing

%% Constants:
A = const.A;
B = const.B;
L = const.L;
density = const.density;

%% Equations:
% Define equations used throughout code in LaTex for use in publishing.

% x_cm:
% x_{cm} = \frac{AB^2 - \alpha\beta^2}{2(AB - \alpha\beta)}

% z_cm:
% z_{cm} = \frac{\frac{1}{2}A^2B - A\alpha\beta + \frac{1}{2}\alpha^2\beta}{AB - \alpha\beta}

% mass:
% m = \rho L(AB -\alpha\beta)

%% Code:
% x_cm
num_x = A*B^2 - alpha*beta^2;
den_x = 2*(A*B - alpha*beta);
x_cm = num_x./den_x;

% z_cm
num_z = A^2*B/2 - A*alpha*beta + alpha^2*beta/2;
den_z = A*B - alpha*beta;
z_cm = num_z./den_z;

% mass
m = density*L*(A*B - alpha*beta);

end