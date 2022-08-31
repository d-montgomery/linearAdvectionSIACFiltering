function [U, t] = oneDLinAdvChebyFilt(x,dt,Tf,N,u0,uB,c,S,nS,nE,cF)
% Solve u_t + c u_x = 0, for x in [-1,1] with filtered initial condition
%       u(x,0) = u0
%       u(boundary,t) = uB

% Inputs:
% x = -cos((0:N)*pi/N) is a column vector for spatial grid
% dt is temporal step size
% Tf is t_final
% a is left endpoint for x 
% b is right endpoint for x
% Nx is number of intervals (Nx+1 gridpoints)
% u0 = @(x) u(x,0) is function with periodic initial condition on (a,b)
% uB = @(t) u(a,t) or u(b,t) is a boundary condition
% c indicates if BC is for left or right boundary 
% S is the filtering matrix
% nS and nE are the starting and ending index for the filter
% cF ensures the ends are accurate (typically = Nd)

% Initial Condition
U = u0(x);

% Filtered IC
Uf = S*U;
% Fix boundary areas
Uf(nS:nS+cF) = U(nS:nS+cF);
Uf(nE-cF:nE) = U(nE-cF:nE);
U = Uf;

% Determine left or right boundary condition
if c < 0 
    n = N+1;
else 
    n = 1;
end


% Get Differentiation Matrix D
D = ChebD(x);

% Initialize t
t = 0;

% TVD Runge-Kutta
while t < Tf
    U(n) = uB(t);
    U1 = U - dt*c*D*(U);
    U1(n) = uB(t+dt);

    U2 = 0.75 * U + 0.25 * (U1 - dt*c*D*U1);
    U2(n) = uB(t+dt/2);

    U = U/3 + 2/3*(U2 - dt*c*D*U2);

    t = t + dt;
end
U(n) = uB(t);