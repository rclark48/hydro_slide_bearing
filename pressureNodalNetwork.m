function [pressure,coord] = pressureNodalNetwork(n,height,const,U,solver)
%% Tittle Block:

% Written by: Reed Clark
% Date Created: 12-31-2018
% Revised By: Reed Clark
% Revision Description: added LaTex equation definitions

%% Description:
% The hydrodynamic pressure is given by the discretized Reynold's
% Equation in 2 dimensions. This function defines the discretization
% and solves for the pressure distribution. The function can solve for
% the nodal pressures in two different ways: one using Gauss-Seidel
% iteration and another just using mldivide and letting MATLAB decide
% how to solve the pressure distribution. Regardless of the solution
% method selected, the function will return the result formatted in
% the same way. The idea of providing two options is to run some test
% cases and see which version of the solver runs faster. The system
% should be diagonally dominant and therefore have a guaranteed
% convergent solution, but the Gauss-Seidel solver may be neccessary
% as the resolution (i.e. the number of nodal pressure equations)
% increases.

%% Inputs:
% n: scalar, fields define the resolution of the nodal network, i.e.
% the step size of the discretization

% height: structure, fields define the inlet height and outlet height

% const: structure, list of shared constants, controlled by the
% constants function

% solver: string, defines which solver to use: Gauss-Seidel, mldivide,
% or both as a sanity check

%% Outputs: 
% pressure: structure.dist, [m], nodal pressure distribution. If the
% mldivide solver is chosen, other fields will be generated within the
% structure as a means to facilitate that solver.

% coord.disc: structure.field, [m], nodal coordinates as disctretized
% by the specified resolutiuon. for code compatness, vectors and
% matrices will be used during computations while the coord array is
% just used to compatly pass information relevant to other functions
% out of the function

%% Constants:
B = const.B;
L = const.L;
eta = const.eta;

hi = height.hi;
ho = heigh.ho;

%% Equations:
% Define equations used throughout code in LaTex for use in publishing.

% Reynold's Equation:
% \frac{\partial }{\partial x} \begin{pmatrix} h^3 \frac{\partial p}{\partial x} \end{pmatrix}  +  \frac{\partial }{\partial y} \begin{pmatrix} h^3 \frac{\partial p}{\partial y} \end{pmatrix} = 6 \eta U \frac{\partial h}{\partial x}

% Simplified Reynold's Equation:
% \begin{pmatrix} 3 h^2 \frac{\mathrm{d} h}{\mathrm{d} x} \end{pmatrix} + h^3 \frac{\partial^2 p}{\partial x^2} + h^3 \frac{\partial^2 p}{\partial y^2} = 6 \eta U \frac{\mathrm{d} h}{\mathrm{d} x}

% dhdy:
% \frac{\partial h}{\partial y} = 0

% Lambda:
% \lambda = \frac{\Delta x}{\Delta y}

% dhdx:
% \frac{\mathrm{d} h}{\mathrm{d} x} = \frac{h_{m,n+1} - h_{m,n-1}}{2 \Delta x}

% dpdx:
% \frac{\partial p}{\partial x} = \frac{p_{m+1,n} - p_{m-1,n}}{2 \Delta x}

% dpdy:
% \frac{\partial p}{\partial y} = \frac{p_{m,n+1} - p_{m,n-1}}{2 \Delta y}

% dp2dx2:
% \frac{\partial^2 p}{\partial x^2} = \frac{p_{m+1,n} + p_{m-1,n} - 2 p_{m,n}}{\Delta x^2}

% dp2dx2:
% \frac{\partial^2 p}{\partial y^2} = \frac{p_{m,n+1} + p_{m,n-1} - 2 p_{m,n}}{\Delta y^2}

% Discretized Reynold's Equation:
% A \begin{pmatrix} p_{m+1,n} - p_{m-1,n} \end{pmatrix} + B \begin{pmatrix} p_{m+1,n} + p_{m-1,n} - p_{m,n} \end{pmatrix} + C \begin{pmatrix} p_{m,n+1} + p_{m,n-1} - p_{m,n} \end{pmatrix} = 6 \eta U \frac{\mathrm{d} h}{\mathrm{d} x}

% A:
% A = \frac{3 h^2}{2 \Delta x} \frac{\mathrm{d} h}{\mathrm{d} x}

% B:
% B = \frac{h^3}{\Delta x^2}

% C:
% C = \frac{h^3}{\Delta y^2}

% Nodal Reynold's Equation (conservation at each interior node):
% C p_{m,n-1} + (B - A) p_{m-1,n} -2 (B + C) p_{m,n} + (A + B) p_{m+1,n} + C p_{m,n+1} = 6 \eta U \frac{\mathrm{d} h}{\mathrm{d} x}
%% Code:
% Grid discretization:
x = linspace(0,B,n); % vector of x-coordinates
coord.disc.x = x;
dx = x(2)-x(1);
m = round(L/B*n); % number of y-direction discretizations
y = linspace(0,L,m); % vector of y-coordinates
coord.disc.y = y;
dy = y(2)-y(1);
lambda = dx/dy; % as currently discretized, lambda should = one
[x_mn, y_mn] = meshgrid(x,y);
coord.disc.x_mn = x_mn;
coord.disc.y_mn = y_mn;

% H
hi_mn = hi.*ones(m,n);
ho_mn = ho.*ones(m,n);
h_mn = hi_mn - (hi_mn - ho_mn)./B.*x_mn;
dhdx = -(hi - ho)./B;

% Network initialization
p_mn = zeros(m,n);

% Solving nodal pressures
switch solver
    case 'Gauss-Seidel'
        pOld = p_mn;
        err = ones(1,2);
        iter = 0;
        
        i = 2:m-1;
        
        j = 1:n;
        jBackward = j(1:end-2);
        jForward = j(3:end);
        jCenter = j(2:end-1);
        
        outerCoeff = 1/(2*(1 + lambda^2));
        plusOne = 1 + 3*(h_mn(i,jForward) - h_mn(i,jBackward))./(4*h_mn(i,jCenter));
        minusOne = 1 - 3*(h_mn(i,jForward) - h_mn(i,jBackward))./(4*h_mn(i,jCenter));
        mixedCoeff = lambda^2;
        speedConst = 3*eta*U*dx*(h_mn(i,jForward) - h_mn(i,jBackward))./(h_mn(i,jCenter).^3);
        
        while all(err(:)) >= const.err.gauss
            % Since dhdy = 0, this special case of applying
            % Gauss-Seidel iteration can actually be vectorially
            % computed.
            plusOneP = plusOne.*p_mn(i,jForward);
            minusOneP = minusOne.*p_mn(i,jBackward);
            mixed = mixedCoeff.*(p_mn(i,jForward) + p_mn(i,jBackward));
            
            p_mn(i,jCenter) = outerCoeff.*(plusOneP + minusOneP + mixed - speedConst);
            
            err = abs((p_mn - pOld)./p_mn);
            
            pOld = p_mn;
            iter = iter + 1;
        end
        pressure.dist = p_mn;
    case 'mldivide'
        i = 1:m; % number of rows
        j = 1:n; % number of columns
        k = m*n; % number of eqations
        rhs = 6*eta*U*dhdx;
        rhs_k = zeros(k,1);
        
        % Spacial pairing of row and column subscripts
        [coord.sub.i, coord.sub.j] = meshgrid(i,j);
        
        % Spacial representation of linear indicies of the nodal
        % network, used to set up the system of equations
        coord.ind = sub2ind([m,n],coord.sub.i,coord.sub.j);
        
        % Build coefficient matrix for pressure system
        pressure.sys = speye(k); % initialize sparse matrix
        for i = 2:m-1
            for j = 2:n-1
                kLeft = coord.ind(i,j-1);
                kUp = coord.ind(i-1,j);
                k = coord.ind(i,j);
                kDown = coord.ind(i+1,j);
                kRight = coord.ind(i,j+1);
                
                A = 3*dhdx*h_mn(k)^2./(2*dx);
                B = h_mn(k)^3./(dx^2);
                C = h_mn(k)^3./(dy^2);
                
                pressure.sys(k,kLeft) = C;
                pressure.sys(k,kUp) = B - A;
                pressure.sys(k,k) = -2*(B + C);
                pressure.sys(k,kDown) = A + B;
                pressure.sys(k,kRight) = C;
                
                rhs_k(k) = rhs;
            end
        end
        p_k = pressure.sys\rhs_k;
        pressure.dist = reshape(p_k,m,n);       
end
end