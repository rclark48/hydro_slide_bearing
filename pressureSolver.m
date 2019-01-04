function [pressure,H] = pressureSolver(U,H,disc,const,solver)
%% Tittle Block:

% Written by: Reed Clark
% Date Created: 12-31-2018
% Revised By: Reed Clark
% Revision Description: renamed function

%% Description:
% The hydrodynamic pressure is given by the discretized Reynold's
% Equation in 2 dimensions. The function can solve for the nodal
% pressures in two different ways: one using Gauss-Seidel iteration
% and another just using mldivide and letting MATLAB decide how to
% solve the pressure distribution. Regardless of the solution method
% selected, the function will return the result formatted in the same
% way. The idea of providing two options is to run some test cases and
% see which version of the solver runs faster. The system should be
% diagonally dominant and therefore have a guaranteed convergent
% solution, but the Gauss-Seidel solver may be neccessary as the
% resolution (i.e. the number of nodal pressure equations) increases.
% Mldivide appears to be the faster solver.

%% Inputs:
% U: scalar, [m/s], bearing speed

% H: structure, fields define the inlet height and outlet height

% disc: structure, [m], nodal coordinates as disctretized by the
% specified resolutiuon. for code compatness, vectors and matrices
% will be used during computations while the coord array is just used
% to compatly pass information relevant to other functions out of the
% function

% const: structure, list of shared constants, controlled by the
% constants function

% solver: string, defines which solver to use: Gauss-Seidel, mldivide,
% or both as a sanity check

%% Outputs: 
% pressure: matrix, [m], nodal pressure distribution. If the
% mldivide solver is chosen, other fields will be generated within the
% structure as a means to facilitate that solver.

% disc: structure, [m], nodal coordinates as disctretized by the
% specified resolutiuon. for code compatness, vectors and matrices
% will be used during computations while the coord array is just used
% to compatly pass information relevant to other functions out of the
% function

%% Constants:
B = const.B;
eta = const.viscosity;

hi = H.hi *1e-6; % [microns --> m]
ho = H.ho *1e-6; % [microns --> m]

n = disc.n;
m = disc.m;
dx = disc.dx;
dy = disc.dy;
lambda = disc.lambda;

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
% H
hi_mn = hi.*ones(m,n);

dhdx = -(hi - ho)./B; % scalar, constant everywhere because of linear profile
h_mn = hi_mn + dhdx.*x_mn;
    H.dhdx = dhdx;
    H.h_mn = h_mn;

% Network initialization
p_mn = zeros(m,n);

% Solving nodal pressures
switch solver
    case 'Gauss-Seidel'
        pOld = p_mn;
        err = ones(1,2);
        iter = 1;
        
        while all(err(:)) > const.err.gauss
            for i = 2:m-1
                for j = 2:n-1
                    outerCoeff = 1./(2*(1 + lambda^2));
                    
                    up = (1 + 3*(h_mn(i,j+1) - h_mn(i,j-1))./(4*h_mn(i,j)))*p_mn(i,j+1);
                    down = (1 - 3*(h_mn(i,j+1) - h_mn(i,j-1))./(4*h_mn(i,j)))*p_mn(i,j-1);
                    leftRight = lambda^2*(p_mn(i+1,j) + p_mn(i-1,j));
                    speed = 3*eta*U*dx*(h_mn(i,j+1) - h_mn(i,j-1))./(h_mn(i,j)^3);
                    
                    p_mn(i,j) = outerCoeff*(up + down + leftRight - speed);
                end
            end
            err = abs((p_mn(2:end-1,2:end-1) - pOld(2:end-1,2:end-1))./p_mn(2:end-1,2:end-1));
            
            pOld = p_mn;
            iter = iter + 1;
        end
        pressure = p_mn;
    case 'mldivide'
        i = 1:m; % number of rows
        j = 1:n; % number of columns
        k = m*n; % number of eqations
        rhs = 6*eta*U*dhdx;
        rhs_k = zeros(k,1);
        
        % Spacial pairing of row and column subscripts
        [disc.sub.i, disc.sub.j] = meshgrid(i,j);
        
        % Spacial representation of linear indicies of the nodal
        % network, used to set up the system of equations
        disc.ind = sub2ind([m,n],disc.sub.i,disc.sub.j);
        
        % Build coefficient matrix for pressure system
        sys = eye(k); % initialize matrix
        for i = 2:m-1
            for j = 2:n-1
                kLeft = disc.ind(i,j-1);
                kUp = disc.ind(i-1,j);
                k = disc.ind(i,j);
                kDown = disc.ind(i+1,j);
                kRight = disc.ind(i,j+1);
                
                A = 3*dhdx*h_mn(k)^2./(2*dx);
                B = h_mn(k)^3./(dx^2);
                C = h_mn(k)^3./(dy^2);
                
                sys(k,kLeft) = C;
                sys(k,kUp) = B - A;
                sys(k,k) = -2*(B + C);
                sys(k,kDown) = A + B;
                sys(k,kRight) = C;
                
                sys = sparse(sys); %convert to sparse for computational efficiency
                
                rhs_k(k) = rhs;
            end
        end
        p_k = sys\rhs_k;
        pressure = reshape(p_k,m,n);       
end
end