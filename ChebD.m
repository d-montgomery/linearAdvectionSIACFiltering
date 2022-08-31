function D = ChebD(x)
    % x must be column vector with Chebyshev nodes
    N = length(x)-1;
    c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
    X = repmat(x,1,N+1);
    dX = X-X';
    D = (c*(1./c)')./(dX +(eye(N+1)));  % off Diagonal Etries
    D = D - diag(sum(D')); %diagonal Entries
end
