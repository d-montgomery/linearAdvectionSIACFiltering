function [U,x,y,t] = twoDLinAdvCheby(dt,Tf,N,u0,uBX,uBY,cx,cy)
% Solve  u_t + cx u_x + cy u_y= 0, for -1 < x,y < 1 and t > 0
%        u(x,y,0) = u0(x,y)
%        u(bndry, y, t) = uBX(y,t)
%        u(x, bndry, t) = uBY(x,t)

% Inputs:
% dt is temporal step size, Tf is t_final
% a is left/lower endpoint for x and y 
% b is right/upper endpoint for x and y
% N is number of intervals (N+1 gridpoints)
% u0 = @(x,y) is u(x,y,0) Initial Condition
% uBX = @(y,t) is u( a or b ,y,t) x Boundary condition
% uBY = @(y,t) is u(x, a or b ,t) y Boundary condition
% cx and cy indicate if BC is for left or right boundary 

% Gauss-Lobatto quadrature nodes on (-1,1)
x = -cos((0:N)*pi/N);
y = x';

U = u0(x,y); % Initial Condition
D = ChebD(x'); % Get Differentiation Matrix D
t = 0; % Initialize t

% Determine left or right boundary condition
if cx < 0  
    nx = N+1;
else
    nx = 1;
end

% Determine lower or upper boundary condition
if cy < 0  
    ny = N+1;
else
    ny = 1;
end

% TVD Runge-Kutta
while t < Tf
    U1 = U - dt*( cx*U*D' + cy*D*U);
    U1(:,nx) = uBX(y,t+dt); % BC updates
    U1(ny,:) = uBY(x,t+dt);

    U2 = 0.75 * U + 0.25 * (U1 - dt*(cx*U1*D' + cy*D*U1));
    U2(:,nx) = uBX(y,t+dt/2); % BC updates
    U2(ny,:) = uBY(x,t+dt/2);

    U = U/3 + 2/3*(U2 - dt*(cx*U2*D' + cy*D*U2));
    
    t = t + dt;
    
    U(:,nx) = uBX(y,t); %BC updates
    U(ny,:) = uBY(x,t);
end