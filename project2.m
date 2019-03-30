


%% need to be passed : r_s'
% Define the confidence level
alpha = 0.95;
% Define the lower and upper bounds to our portfolio
lb = [-inf(n,1); zeros(S,1); -inf ];
ub = [];
% Define the inequality constraint matrices A and b
A = [ -r_s' -eye(S) -ones(S,1) ];
b = [ zeros(S, 1)];
% Define the equality constraint matrices A_eq and b_eq
Aeq = [ ones(1,n) zeros(1,S) 0 ];
beq = 1;

% Define our objective linear cost function c
k = (1 / ( (1 - alpha) * S) );
c = [ zeros(n,1); k * ones(S,1); 1 ; -0.1*mu' zeros(1,S)];
% Set the linprog options to increase the solver tolerance
options = optimoptions('linprog','TolFun',1e-9);

% Use 'linprog' to find the optimal portfolio
y = linprog( c, A, b, Aeq, beq, lb, ub, [], options );

% Retrieve the optimal portfolio weights
x = y(1:n);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 2: Ellipsoidal uncertainty set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Uncertainty set size
Theta = diag( diag(Q) ) ./ N;

% Square root of Theta
sqrtTh = sqrt(Theta);

% Confidence level
alpha = 0.9;

% Scaling parameter epsilon for uncertainty set
ep = sqrt( chi2inv(alpha, n) );

%% Setup our input parameters for fmincon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   min     gamma + (1 / [(1 - alpha) * S]) * sum( z_s ) - lamda * mu'x
%   s.t.    z_s   >= 0,                 for s = 1, ..., S
%           z_s   >= -r_s' x - gamma,   for s = 1, ..., S  
%           1' x  =  1,







%          