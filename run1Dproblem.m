% Solve 1D Advection -----------------------------------------------------
% u_t + c u_x = 0, -1 < x < 1, t > 0
% u(x,0) = 0.2 sin(10x) - 0.5  if x <= -0.25
%        = 0.2 sin(10x) + 0.5  if x > -0.25
% u(-1,t) = 0.2 sin(10*(-1-t)) - 0.5 
%-------------------------------------------------------------------------
clear 
close all

% --- Flags --------------------------------------------------------------
% saveWS = 1   : saves the unfiltered data U, x and N
% saveWS = 2   : saves the filters

% loadFilt = 1 : loads S matrix with  m = 3; k = 8; Nd = 13;
%                This will drastically speed up test = 2

% test = 1     : get unfiltered data via Chebyshev Collocation
% test = 2     : get filtered solution

saveWS = 0; % this will save in main folder, not data
loadFilt = 1; % loads from data folder (adjust as needed)
test = 2;

% Set up N values for spatial grid
cs = 0; % 0 
cf = 4; % 4 gives N = [100, 150, 200, 250, 300]
j = cf - cs; % for initializing errors and loops
N = 100 + 50*(cs:cf);

% Temporal Parameters
Tf = 1;
dt = 1e-5;

% Initial Condition
u0 = @(x) 0.2*sin(10*(x-1))-0.5 +...
          heaviside(x+0.25).*(0.2*sin(10*(x-1)) + ...
          0.5 - (0.2*sin(10*(x-1)) - 0.5));    

% Boundary Condition
uB = @(t) 0.2*sin(10*((-1-t)-1)) - 0.5;
c = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 1 - Get Unfiltered Data via Chebyshev Collocation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if test == 1
    
    % Initialize Error and cell to store data for later
    Err = zeros(j+1, 1);
    UXN = cell(j+1,4); %UXN = {U, x, N} for each N value
    AbsErrUnFilt = cell(j+1,4);

    % Open a Figure for Solutions
    color = ['b'; 'g'; 'r'; 'c'; 'm'];
    fig1 = figure(1);
    axis([-1,1, -1,1])
    p=gca;

    
    for k = 1:j+1

        % Get Approxomition
        [U, x, tf] = oneDLinAdvCheby(dt, Tf, N(k), u0, uB, c);

        ue = u0(x-c*tf);
        AbsErr = abs(ue - U);
        AbsErrUnFilt{k} = AbsErr;
        Err(k) = norm(AbsErr,inf)/norm(ue,inf);

        % Store U, x, N, and Tf for filtering later
        UXN{k,1} = U; UXN{k,2} = x; UXN{k,3} = N(k); UXN{k,4} = Tf;

        % Plot the Error
        subplot(1,2,2)
        semilogy(x, AbsErr, color(k),'DisplayName', ...
            ['N = ',num2str(N(k))])
        xlabel('$x$', 'Interpreter', 'Latex', 'FontSize', 18)
        ylabel('$error$', 'Interpreter', 'Latex', 'FontSize', 18)
        hold on

        % Plot the Functions
        if k == 1
            xf = -cos((0:N(end))'*pi/N(end));
            subplot(1,2,1)
            plot(xf,u0(xf-c*tf),'-.', 'MarkerSize', 4, 'DisplayName', 'Exact')
            hold on
        end

        subplot(1,2,1)
        plot(x,U, color(k), 'DisplayName', ...
            ['N = ',num2str(N(k))])
        axis([-1,1, -1,1])
        xlabel('$x$', 'Interpreter', 'Latex', 'FontSize', 18)
        ylabel('$u$', 'Interpreter', 'Latex', 'FontSize', 18)
        hold on

    end
    hold off
    subplot(1,2,1)
    legend('Location', 'NorthWest', 'FontSize', 14)
    
    subplot(1,2,2)
    legend('Location', 'NorthWest', 'FontSize', 14)

    % Save U, x, N, and Tf for filtering later
    if saveWS == 1
            save('UXN_1D.mat', 'UXN')
            save('ErrUnfiltered1D.mat','AbsErrUnFilt')
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 2 - Filter Inititial Solution at each time step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if test == 2
    
    % Parameters for filter
    m = 3; k = 8;
    Nd = 13;
    cF = Nd;
    filters = cell(j+1,3);
    if loadFilt == 1 % get pre-computed filters with m = 3; k = 8; Nd = 13
        load('data/filters.mat')
    end

    % Load Unfiltered data for comparison
    load('data/UXN_1D.mat')

    % Store Errors
    Err = cell(j+1, 1);
    
    % Open a Figure for Solutions
    color = ['b'; 'g'; 'r'; 'c'; 'm'];
    fig2 = figure(2);
    axis([-1,1, -1,1])
    
    % Store Filtered Data
    UfXN = cell(j+1,3);

    ctr = cs;
    for i = 1:j+1
        ctr = ctr + 1;
      
        % Gauss-Lobatto quadrature nodes on (-1,1)
        x = -cos(pi*(0:N(i))/N(i))'; % creates x from - to +
        
        % Get filter 
        if ~loadFilt
            [S, nS, nE] = filtering.matrix(x, m, k, Nd);
            cF = ceil(Nd);
            if saveWS == 1
                filters(i,:) = {S,nS,nE};
            end
        else
            S = filters{ctr,1}; nS = filters{ctr,2}; nE = filters{ctr,3};
        end
        
        % Get Unfiltered Approxomition
        U = UXN{ctr,1};
    
        % Get Filtered Approximation
        [Uf, tf] = oneDLinAdvChebyFilt(x,dt,Tf,N(i),u0,uB,c,S,nS,nE,cF);
        UfXN{i,1} = Uf; UfXN{i,2} = x; UfXN{i,3} = N(i);
        
        ue = u0(x-c*tf);
        Err{i} = abs(ue - Uf);
        
        subplot(1,2,1)
        plot(x,U, color(ctr), 'DisplayName', ...
            ['N = ',num2str(N(i))])
        axis([-1,1, -1,1])
        xlabel('$x$', 'Interpreter', 'Latex', 'FontSize', 18)
        ylabel('$u$', 'Interpreter', 'Latex', 'FontSize', 18)
        title('Unfiltered Data', 'FontSize', 18)
        hold on
        
        subplot(1,2,2)
        plot(x,Uf, color(ctr), 'DisplayName', ...
            ['N = ',num2str(N(i))])
        axis([-1,1, -1,1])
        xlabel('$x$', 'Interpreter', 'Latex', 'FontSize', 18)
        ylabel('$u$', 'Interpreter', 'Latex', 'FontSize', 18)
        title('Filtered Data','FontSize', 18)
        hold on        
    end
    subplot(1,2,1)
    legend('Location', 'NorthWest', 'FontSize', 16)
    subplot(1,2,2)
    legend('Location', 'NorthWest', 'FontSize', 16)
    hold off
    
    % Side-by-Side Error Plots
    figure(3) 
    load('data/ErrUnfiltered1D.mat')
    color = ['b'; 'g'; 'r'; 'c'; 'm'];
    p=gca;
    ctr = cs;
    for i = 1:j+1
        ctr = ctr+1;
        x = UfXN{i,2}; Uf = UfXN{i,1};
        errU = AbsErrUnFilt{ctr};
        
        subplot(1,2,1)
        semilogy(x,errU, color(ctr), 'DisplayName', ...
            ['N = ',num2str(N(i))])
        xlabel('$x$', 'Interpreter', 'Latex', 'FontSize', 18)
        ylabel('$error$', 'Interpreter', 'Latex', 'FontSize', 18)
        title('Unfiltered Error','FontSize', 18)
        hold on 
        
        subplot(1,2,2)
        semilogy(x,Err{i}, color(ctr), 'DisplayName', ...
            ['N = ',num2str(N(i))])
        xlabel('$x$', 'Interpreter', 'Latex', 'FontSize', 18)
        ylabel('$error$', 'Interpreter', 'Latex', 'FontSize', 18)
        title('Filtered Error','FontSize', 18)
        hold on 
    end
    subplot(1,2,1)
    legend('Location', 'NorthWest', 'FontSize', 16)
    
    % Save Filters for Later
    if saveWS == 1
            save('filters.mat', 'filters')
    end
end
