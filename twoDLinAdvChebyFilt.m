function [U] = twoDLinAdvChebyFilt(x,dt,Tf,N,u0,uXB,uYB,cx,cy,S,nS,nE,cF)
% Solve with filtered initial condition 
%        u_t + cx u_x + cy u_y= 0, for -1 < x,y < 1 and t > 0 
%        u(x,y,0) = u0(x,y)
%        u(bndry, y, t) = uBX(y,t)
%        u(x, bndry, t) = uBY(x,t)
% x = -cos((0:N)*pi/N) is a row vector for spatial
% dt is temporal step size, Tf is t_final
% N is number of intervals (N+1 gridpoints)
% u0 = @(x,y) is u(x,y,0) Initial Condition
% uXB = @(y,t) is u( bndry ,y,t) x Boundary condition
% uYB = @(y,t) is u(x, bndry ,t) y Boundary condition
% cx and cy indicate if BC is for left or right boundary 
% S is the filtering matrix
% nS and nE are the starting and ending index for the filter
% cF ensures the ends are accurate (typically = Nd)

% Gauss-Lobatto quadrature nodes on (-1,1)
y = x';

D = ChebD(x'); % Get Differentiation Matrix D
t = 0; % Initialize t
U = u0(x,y); % Initial Condition
Uf = S*U*S'; % Filtered Initial Condition

% Ensure area around boundaries is accurate
Uf(nS , nS:nS+cF) = U(nS , nS:nS+cF);
Uf(nS:nS+cF , nS ) = U(nS:nS+cF , nS);
Uf(nE, nE-cF:nE) = U(nE, nE-cF:nE);
Uf(nE-cF:nE,nE) = U(nE-cF:nE,nE);
U = Uf;

% Determine left or right boundary condition
if cx < 0  
    nx = N+1;
else
    nx = 1;
end
% Determine bottom or top boundary condition
if cy < 0  
    ny = N+1;
else
    ny = 1;
end

% TVD Runge-Kutta
while t < Tf
    U1 = U - dt*( cx*U*D' + cy*D*U);
    U1(:,nx) = uXB(y,t+dt); % BC updates
    U1(ny,:) = uYB(x,t+dt);

    U2 = 0.75 * U + 0.25 * (U1 - dt*(cx*U1*D' + cy*D*U1));
    U2(:,nx) = uXB(y,t+dt/2); % BC updates
    U2(ny,:) = uYB(x,t+dt/2);

    U = U/3 + 2/3*(U2 - dt*(cx*U2*D' + cy*D*U2));
    
    t = t + dt;
    U(:,nx) = uXB(y,t); %BC updates
    U(ny,:) = uYB(x,t);
end