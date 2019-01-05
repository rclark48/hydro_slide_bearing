function [force, moment] = integrator(pressure,U,disc,H,const,com)
%% Tittle Block:

% Written by: Reed Clark
% Date Created: 1-4-2019
% Revised By: Reed Clark
% Revision Description: corrected year in date created

%% Description:
% This function numerically integrates the pressure distribution to
% obtain the hydrodynamic load support forces and moments.

%% Inputs:
% pressure: matrix, [m], nodal pressure distribution

% U: scalar, [m/s], bearing trasnlation speed

% disc: structure, [m], nodal coordinates as disctretized by the
% specified resolutiuon. for code compatness, vectors and matrices
% will be used during computations while the coord array is just used
% to compatly pass information relevant to other functions out of the
% function

% H: structure, fields define the inlet height and outlet height

% const: structure, list of shared constants, controlled by the
% constants function

% com: structure, contains center of mass coordinates

%% Outputs: 
% force.loadSupport: scalar, [N], value of hydrodynamic load support

% force.drag: scalar, [N], hydrodynamic drag

% moment.loadSupport: scalar, [Nm], moment about com due to
% hydrodynamic load support.

%% Constants:
eta = const.viscosity;

dhdx = H.dhdx;
h_mn = H.h_mn;

x = disc.x;
m = disc.m;
x_mn = disc.x_mn;
y = disc.y;

x_cm = com.x;

%% Equations:
% Define equations used throughout code in LaTex for use in publishing.

% Hydrodynamic load support:
% \int_{0}^{L} \int_{0}^{B} p(x,y) dx dy

% Hydrodynamic drag:
% -\frac{1}{2} \frac{\mathrm{d} h}{\mathrm{d} x} \int_{0}^{L} \int_{0}^{B} p(x,y) dx dy - \eta U \int_{0}^{L} \int_{0}^{B} \frac{1}{h(x)} dx dy

% Hydrodynamic load support moment:
% M = \int_{0}^{L} \{ \int_{0}^{x_{cm}} p(x,y)x dx \int_{x_{cm}}^{B} p(x,y)x dx \} dy

% Hydrodynamic load support moment (discretized version):
% M = \int_{0}^{L} \int_{0}^{x_{cm}} p(x,y) x dx dy - \int_{0}^{L} \int_{x_{cm}}^{B} p(x,y) x dx dy

%% Code:
% force.loadSupport:
F = trapz(x,trapz(y,pressure));
force.loadSupport = F;

% force.drag:
D = -dhdx./2*F - eta*U*trapz(x,trapz(y,1./h_mn));
force.drag = D;

% moment.loadSupport:
ind = find(x > x_cm);
xL = [x_mn(:,1:ind(1) - 1) x_cm*ones(m,1)]; % nodes <= x_cm
xU = [x_cm*ones(m,1) x_mn(:,ind(1):end)]; % nodes >= x_cm

p_cm = interp1(x,pressure',x_cm)';

pL = [pressure(:,1:ind(1) - 1) p_cm]; % pressure dist below x_cm
pU = [p_cm pressure(:,ind(1):end)]; % pressure dist above x_cm

mPlus = trapz(xL(1,:),trapz(y,pL.*xL));
mMinus = trapz(xU(1,:),trapz(y,pU.*xU));
M= mPlus - mMinus;
moment.loadSupport = M;
end