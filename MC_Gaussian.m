function [S] = MC_Gaussian(initPrice,mu_FF, Q_FF)

% %%%%%%%%%%TEST PARAMETERS DELETE THIS CODE ONCE FINISHED TESTING
% % Stock  parameters (weekly)
% initPrice    = [100; 150];          % Initial price of stock A
% mu_FF    = [0.1 / 100; 0.2/100];    % expected return
% Q_FF = [0.000064 0.0000576; 0.0000576 0.000144];    % covariance matrix
% 
% 
% %%%%%%%%%%TEST PARAMETERS DELETE THIS CODE ONCE FINISHED TESTING

%This function calculates stock price scenarios using the Gaussian Monte
%Carlo simulations

%Calculating correlation matrix
corrMat = corrcov(Q_FF);

% Experimental parameters
numbYear = 3; %Number of year that will be simulated
T   = 52 * numbYear;       % Time window  
N   = 52 * numbYear;      % Number of steps (half a week per time step)
dt  = T / N;    % Timestep 
n =  size(mu_FF,1);%number of assets
sigmaQ = sqrt(diag(Q_FF)); %standard deviations of assets

% Number of simulated price paths
nPaths = 5000;

%--------------------------------------------------------------------------
% 2.3 Discrete-time approach using the Matlab function 'mvnrnd()'
% This process behaves in the exact same way as the one above.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find our correlation transformation matrix L
L = chol(corrMat, 'lower');

% Allocate space for our simulations for each asset
S = cell(n,1);

%S_A(1,:) = S0_A;   
%S_B(1,:) = S0_B; 


for k = 1: n;
    %Allocate space for each asset in simulation
   S_temp = zeros(N+1, nPaths); 
   % Set initial price for assets
   S_temp(1,:) = initPrice(k);
   S{k} = S_temp;
end
S{k}(2,1)
for i = 1:nPaths
    for j = 1:N
        
        % Convert the independent RVs to correlated RVs
        xi = L * randn(2, 1); 
       
        for k = 1 : n
        
        %S{k}(j+1, i) = S{k}(j, i) * exp( ( mu_FF(k) - 0.5 * (sigmaQ(k))^2 ) * dt ...
        %                + sigmaQ(k) * sqrt(dt) * xi(k) );
        S{k}(j+1, i) = S{k}(j, i) + S{k}(j, i) * ...
                                        mvnrnd(dt * mu_FF(k), dt * (sigmaQ(k))^2);
        
        %discS_A2(j+1, i) = discS_A2(j, i) + discS_A2(j, i) * ...
         %                               mvnrnd(dt * mu_A, dt * sigma_A^2);
        end  
    end
end 


%--------------------------------------------------------------------------
% 2.4 Plot the paths of all the simulations 
%--------------------------------------------------------------------------
fig2 = figure(2);
plot(0:N, S{1})
title('Stock A Price Evolution (Discrete)', 'FontSize', 14)
ylabel('Stock Price','interpreter','latex','FontSize',12);
xlabel('Time','interpreter','latex','FontSize',12);
xlim([0 N])