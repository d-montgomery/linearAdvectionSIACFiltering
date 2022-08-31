% Class that contains all functions necessary for filter
classdef filtering
methods(Static)

% --- Pmk Polynomials-----------------------------------------------------
function [Pmk] = polynomialsPmk(m,k, handleFlag)
    % Compute the Pmk polynomials symbolically
    % If handleFlag = 1, the function returns Pmk = @(x) as fuction handle
    syms x

    % Weight function
    wk = (1 - x^2)^(k+1);

    % Weighted Norm and Inner Product
    normWk = @(f) sqrt(int(f^2*wk, x, -1,1));
    innerWk = @(f,g) int(f*g*wk, x, -1,1);

    % Generate {r_j} orthonormal basis of span({x^j}), j = 1, ..., m 
    % using the Gram-Schmidt orthogonalization procedure.
    r = cell(m,1);
    r{1} = x/normWk(x);
    for j = 1:m-1
        S = 0;
        for i = 1:j
            S = S + r{i}*innerWk(x^(j+1), r{i});
        end
        v = x^(j+1) - S;
        r{j+1} = simplify(v/normWk(v));
    end

    % Compute Q(x)
    Q = 0;
    for j = 1:floor(m/2)
        Q = Q + r{2*j} * innerWk(1,r{2*j});
    end
    Q = 1 - Q;

    % Compute c
    c = 1/innerWk(Q,1);

    % Get Polynomial P^{m,k}(x)
    Pmk = simplify(c*wk*Q);
    if handleFlag == 1
        Pmk = matlabFunction(Pmk); % returns Pmk = @(x) ...
    end
end % pmk polynomial function

       


% --- Clenshaw-Curtis Quadrature -----------------------------------------
function I = clenshaw_curtis(f,Q,xj,es) 
    % f = @(t) .... function to be integrated (polynomial of degree p)
    % Q is number of Chebychev Points (choose Q + 1 = p,)
    % c is center of integral 
    % es is distance left and right from center.  
    %         i.e. a = xj - es, b = xj + es

    t = xj-es*cos(pi*(0:Q)'/Q);  
    fx = es*f(t)/(2*Q) ; 
    g = real(fft(fx([1:Q+1, Q:-1:2]))); 
    a = [g(1) ; g(2:Q)+g(2*Q:-1:Q+2) ; g(Q+1)] ; 
    w = 0*a'; w(1:2:end) = 2./(1-(0:2:Q).^2); 
    I = w*a;
end % clenshaw_curtis


% --- Lagrange Polynomials -----------------------------------------------
function l = LagrangeL(i, x)
    % i -> l_i
    % x -> data
    syms t
    [m,n] = size(x);  
    N = max(m,n);
    l = 1;
    for k = 1:N
        if k ~= i
           l = l.*(t - x(k))/(x(i) - x(k));  
        end
    end   
    l = matlabFunction(l);
end
       

% Filtering matrix
function [S, nLeft, nRight] = matrix(x, m, k, Nd)
    % Inputs
    % x is the discretized domain used in numerical scheme
    % m and k are the associated Dirac-Delta parameters
    % Nd is number of Chebyshev points spanned by the kernel

    % Outputs
    % S = filter matrix
    % nLeft = index where filter starts
    % nRight = index where filter stops

    N = length(x) - 1;

    % Number of quadrature nodes in compact support domain
    Q = 2*(floor(m/2) + k+1) + N+1;

    % Initialize S
    S = eye(N+1,N+1);

    % Calculate epsilon es
    es = sin( pi*Nd / (2*N) );

    % Calculate Dirac-Delta polynomial
    Pmk = filtering.polynomialsPmk(m,k,1);
    delta = @(t) heaviside(t+es).*Pmk(t./es)./es - ...
        heaviside(t-es).*Pmk(t./es)./es;

    % Determine where filter starts and stops 
    % (i.e. don't filter near boundaries)
    nLeft = ceil(N/pi * acos(1 - es) );
    nRight = floor(N/pi * acos(es - 1));

    % Calculate each entry of S
    for i = nLeft:nRight
        l = filtering.LagrangeL(i,x);
        for j = nLeft:nRight
            f = @(t) l(t).*delta(x(j)-t);
            S(i,j) = filtering.clenshaw_curtis(f, Q, x(j), es);
        end
    end
    S(nLeft,nLeft) = 1;
    S(nRight,nRight) = 1;
end
       

end %Static
end %class