% Solve 2D Advection -----------------------------------------------------
% u_t + cx u_x + cy u_y = 0, for -1 < x,y < 1 , t >0
% u(x,y,0) = cos( 4pi*sqrt( x^2 + (y+0.5)^2 ) ) if x^2 + (y+0.5)^2 <= 1/16
%          = 0  otherwise 
%-------------------------------------------------------------------------
clear 
close all

% --- Flags -------------------------------------------------
% saveWS = 1   : saves the unfiltered data U, x and N
% saveWS = 2   : saves the filters

% loadFilt = 1 : loads S matrix with  m = 3; k = 8; Nd = 13;
%                This will drastically speed up test = 2

% test = 1     : get unfiltered data via Chebyshev Collocation
% test = 2     : get filtered solution

% surfFlag = 1 : surf of exact soln and unfiltered soln (N = 150 is best)
% surfFlag = 2 : surf of unfiltered soln and filtered soln (Uses N = 150)

saveWS = 0; % this will save in main folder not data
loadFilt = 1; % loads from data folder (adjust as needed)
test = 2;
surfFlag = 0;

% Set up N values for spatial grid
cs = 0; % 0 
cf = 1; % 4 gives N = [100, 150, 200, 250, 300]
j = cf - cs; % for initializing errors and loops
N = 100 + 50*(cs:cf);

% Temporal Parameters
Tf = .5;
dt = 1e-5;

% Initial Condition
u0 = @(x,y) cos(4*pi*sqrt(x.^2 + (y+0.5).^2)) +...
    - heaviside(x.^2+(y+0.5).^2-1/16).*cos(4*pi*sqrt(x.^2+(y+0.5).^2));    

% Boundary Conditions
uBX = @(y,t) u0(-1,y-t);
cx = 1;

uBY = @(x,t) u0(x-t,-1);
cy = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 1 - Get Data via Chebyshev Collocation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if test == 1
    
    % Initialize Error and cell to store data for later
    Err = zeros(j+1, 1);
    UXN = cell(j+1,6); %UXN = {U, x, N, tf, UmidY, AbsErr} for each N value

    % Open a Figure for Solutions
    color = ['b'; 'g'; 'r'; 'c'; 'm'];
    fig1 = figure(1);
    axis([-1,1, -1.1,1.1])
    p=gca;

    for k = 1:j+1

        % Get Approxomition
        [U,x,y,tf] = twoDLinAdvCheby(dt,Tf,N(k),u0,uBX,uBY,cx,cy);
        UmidY = U(1+N(k)/2,:); % Slice at y = 0
        
        % Calculate Error
        ue = u0(x-cx*tf, y(1+N(k)/2) -cy*tf);
        AbsErr = abs(ue - UmidY);
        Err(k) = max(AbsErr)/max(ue);
        
        % Store U and x for filtering later
         UXN{k,1} = U; UXN{k,2} = x; UXN{k,3} = N(k); UXN{k,4} = tf; 
         UXN{k,5} = UmidY; UXN{k,6} = AbsErr;

        % Plot the Error
        subplot(1,2,2)
        semilogy(x, AbsErr, color(k))
        xlabel('$x$', 'Interpreter', 'Latex', 'FontSize', 16)
        ylabel('$error$', 'Interpreter', 'Latex', 'FontSize', 16)
        hold on

        % Plot the Exact Solution for largest N
        if k == 1
            xf = -cos((0:N(end))*pi/N(end));
            subplot(1,2,1)
            plot(xf, u0(xf-cx*tf, xf(1+N(end)/2)-cy*tf), ...
                'DisplayName', 'Exact')
            hold on
        end
        
        % Plot the Approximations 
        subplot(1,2,1)
        plot(x, UmidY,color(k), 'DisplayName',...
            ['N = ',num2str(N(k))])
        axis([-1,1, -1.1,1.1])
        xlabel('$x$', 'Interpreter', 'Latex', 'FontSize', 16)
        ylabel('$u$', 'Interpreter', 'Latex', 'FontSize', 16)
        hold on

    end
    hold off
    legend('Location', 'NorthWest')

    % Store U and x for filtering later
        if saveWS == 1
            fileName = 'UXN_2D.mat';
                save(fileName, 'UXN');
        end
    
    % Surf of Exact Solution and Numerical Solution
    if surfFlag == 1
        figure(2)
        c = hsv; % nice color scheme
        
        % Surf of Exact Solution
        subplot(1,2,1)
        surf(x,y,u0(x-cx*tf, y -cy*tf))
        colormap(c)
        caxis([-1 1])
        axis([-1,1,-1,1,-1,1])
        xlabel('x','FontSize',16)
        ylabel('y','FontSize',16)
        zlabel('u','FontSize',16)
        title(['Exact Solution at t = ', num2str(tf)], 'FontSize',18)
    
        % Surf of Approximation
        subplot(1,2,2)
        surf(x,y,U)
        axis([-1,1,-1,1,-1,1])
        colormap(c)
        caxis([-1 1])
        xlabel('x','FontSize',16)
        ylabel('y','FontSize',16)
        zlabel('u','FontSize',16)
        title(['Numerical Solution at t = ', num2str(tf),...
            ' with N = ', num2str(N(end))],'FontSize',18)
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 2 - Filter All Solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if test == 2
    
    % Parameters for Filter
    m = 3; k = 8;
    Nd = 13;
    cF = Nd;
    filters = cell(j+1,3);
    if loadFilt == 1 % get pre-computed filters with m = 3; k = 8; Nd = 13
        load('data/filters.mat')
    end
    
    % Load Unfiltered data for comparison
    load('data/UXN_2D.mat')

    % Store Filtered Data
    UFXN = cell(j+1,3);
    
    % Open a Figure for Solutions
    color = ['b'; 'g'; 'r'; 'c'; 'm'];
    fig1 = figure(1);
    axis([-1,1, -1,1])
    
    % Counter from choich of N value
    ctr = cs;

    for i = 1:j+1
        ctr = ctr +1;
        
        % Unfiltered data
        U = UXN{ctr,1}; 
        x = UXN{ctr,2};
        y = x';
        Tf = UXN{ctr,4};
        UmidY = UXN{ctr,5};
        errU = UXN{ctr,6};
        
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
        
        % Get Filtered Approx
        Uf =twoDLinAdvChebyFilt(x,dt,Tf,N(i),u0,uBX,uBY,cx,cy,S,nS,nE,cF);
        UfMidY = Uf(1+N(i)/2,:); % solution for fixed y = 0
        
        % Calculate Error
        ue = u0(x-cx*Tf, y(1+N(i)/2) -cy*Tf);
        AbsErr = abs(ue - UfMidY);
        
        % Store Filted Approx (for surf plot)
        UFXN{ctr,1} = Uf; 
        UFXN{ctr,2} = UfMidY;
        UFXN{ctr,3} = AbsErr;
        
        
        subplot(1,2,1)
        plot(x,UmidY, color(ctr), 'DisplayName', ...
            ['N = ',num2str(N(i))])
        axis([-1,1, -1,1])
        xlabel('$x$', 'Interpreter', 'Latex', 'FontSize', 18)
        ylabel('$u$', 'Interpreter', 'Latex', 'FontSize', 18)
        title('Unfiltered Data','FontSize', 18)
        hold on
        
        subplot(1,2,2)
        plot(x,UfMidY, color(ctr), 'DisplayName', ...
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
    
    % Store Uf and x for later
    if saveWS == 1
        fileName = 'UFXN_2D.mat';
            save(fileName, 'UFXN');
    end
    
    % Open a Figure for Error Plot
    color = ['b'; 'g'; 'r'; 'c'; 'm'];
    fig2 = figure(2);
    axis([-1,1, -1,1])
    p=gca;
    ctr = cs;
    for i = 1:j+1
        ctr = ctr+1;
        x = UXN{ctr,2}; errUf = UFXN{ctr,3};
        errU = UXN{ctr,6};
        
        subplot(1,2,1)
        semilogy(x,errU, color(ctr), 'DisplayName', ...
            ['N = ',num2str(N(i))])
        xlabel('$x$', 'Interpreter', 'Latex', 'FontSize', 18)
        ylabel('$error$', 'Interpreter', 'Latex', 'FontSize', 18)
        title('Unfiltered Error','FontSize', 18)
        hold on 
        
        subplot(1,2,2)
        semilogy(x,errUf, color(ctr), 'DisplayName', ...
            ['N = ',num2str(N(i))])
        xlabel('$x$', 'Interpreter', 'Latex', 'FontSize', 18)
        ylabel('$error$', 'Interpreter', 'Latex', 'FontSize', 18)
        title('Filtered Error','FontSize', 18)
        hold on 
    end
    subplot(1,2,1)
    legend('Location', 'NorthWest', 'FontSize', 16)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Surf of Unfiltered and Filtered Solutions (test ~= 1 or 2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if surfFlag == 2
    load('data/UXN_2D.mat')
    load('data/UXNfilt_2D.mat')

    x = UXN{2,2};
    U = UXN{2,1};
    Uf = UFXN{2,1};
    N = UXN{2,3};

    figure(2)
    c = hsv;
    subplot(1,2,1)
    surf(x,x,U)
    axis([-1,1, -1,1, -1,1])
    colormap(c)
    caxis([-1 1])
    xlabel('x', 'FontSize',16)
    ylabel('y', 'FontSize',16)
    zlabel('u', 'FontSize',16)
    title(['Unfiltered Solution at t = ', num2str(Tf),...
        ' with N = ', num2str(N)], 'FontSize',18)

    subplot(1,2,2)
    surf(x,x,Uf)
    axis([-1,1, -1,1, -1,1])
    colormap(c)
    caxis([-1 1])
    xlabel('x', 'FontSize',16)
    ylabel('y', 'FontSize',16)
    zlabel('u', 'FontSize',16)
    title(['Filtered Solution at t = ', num2str(Tf),...
        ' with N = ', num2str(N)], 'FontSize',18)
end