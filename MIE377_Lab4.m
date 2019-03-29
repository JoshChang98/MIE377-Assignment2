%% MIE377 (Winter 2019) - Laboratory 4
% The purpose of this program is to solve a robut optimization problem to 
% construct robust portfolios with an ellipsoidal uncertainty structure. We 
% will use the MATLAB "fmincon" function to solve this problem.
%
% TA: Ricardo Pillaca
% Instructor: Giorgio Costa

clc
clear all
format short

% Program Start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 1: Data pre-processing 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the sample historical data
load('lab2data.mat')

% Calculate the asset and factor returns (factor models use returns, not
% prices)
rets    = prices( 2:end, : ) ./ prices( 1:end - 1, : ) - 1;
facRets = sp500price( 2:end , 1 ) ./ sp500price( 1:end - 1, 1 ) - 1;

% Number of assets
n = size(rets,2);

% Number of observations;
N = size(rets, 1);

% Calculate the factor expected excess return from historical data using
% the geometric mean
mu = (geomean(rets + 1) - 1)';

% Calculate the asset covariance matrix
Q = cov(rets);

% Risk aversion parameter
lambda = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 2: Ellipsoidal uncertainty set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Uncertainty set size
Theta = diag(diag(Q))./N; 

% Square root of Theta
sqrtTh = sqrt(Theta);

% Confidence level
alpha = 0.9;

% Scaling parameter epsilon for uncertainty set
ep = sqrt(chi2inv(alpha,n));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 3: Setup our input parameters for fmincon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% min   lambda * (x' * Q x) - mu' x + epsilon * norm (sqrtTh * x)
% s.t.  sum(x) == 1
%       x >= 0

%-------------------------------------------------------------------------- 
% 3.1 Inequality constraints:
% fmincon accepts linear inequality constraints of the form "A x <= b".
%--------------------------------------------------------------------------

% The only inequality constraint is the one on short selling. However, this
% can be applied as a bound on our variable x.

% Linear inequality Constraint bounds
b = [];
A = []; 

% Lower and upper bounds on variables
lb = zeros(n,1);
ub = ones(n,1);

%--------------------------------------------------------------------------
% 3.2 Equality constraints: 
% We only have the budget constraint
%--------------------------------------------------------------------------

beq = 1;
Aeq = ones(1,n);

%--------------------------------------------------------------------------
% 3.3 Initial solution: 
% fmincon requires an initial feasible solution 
%--------------------------------------------------------------------------
% Define an initial portfolio ("equally weighted" or "1/n portfolio")
x0 = repmat(1/n,n,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 4: For comparison, find the nominal portfolio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Allocate space for both our nominal and robust solutions
x = zeros(n, 2);

% Increase the tolerance of 'quadprog'
options = optimoptions('quadprog','TolFun',1e-9);

% Find the nominal portfolio with 'quadprog'
x(:,1) = quadprog( lambda * 2 * Q, -mu, A, b, Aeq, beq, lb, ub, [], options );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 5: Start 'fmincon' solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solve using fmincon to find our optimal portfolio weights 
x(:,2) = fmincon(@(x)objFun(x, mu, Q, lambda, sqrtTh, ep), x0, A, b, ...
                Aeq, beq, lb, ub, @(x)nonlcon(x));
        
% Print both solutions to the console
display(x)

% Plot our portfolio weights for comparison
fig1 = figure(1);
bar(x','stacked');
set(gca,'TickLabelInterpreter', 'latex','fontsize',22);
set(gca, 'XTickLabel', {'Nomimal','Robust'},'fontsize',22);
ylabel('Weight','interpreter', 'latex','FontSize',24);
title('Portfolio Weights','interpreter', 'latex','FontSize',24);

set(fig1,'Units','Inches', 'Position', [0 0 10, 10]);
    pos3 = get(fig1,'Position');
    set(fig1,'PaperPositionMode','Auto','PaperUnits','Inches',...
        'PaperSize',[pos3(3), pos3(4)])

%-------------------------------------------------------------------------- 
% 5.1 Define the objective function:
% We must specify our nonlinear objective as a separate function. fmincon
% accepts the "norm( )" function.
%--------------------------------------------------------------------------

function f = objFun(x, mu, Q, lambda, sqrtTh, ep)

     f = lambda * (x' * Q * x) - mu'* x + ep * norm (sqrtTh * x);

end

%-------------------------------------------------------------------------- 
% 5.2 Define the equality and inequality nonlinear constraints:
% fmincon accepts nonlinear constraints, but these must be defined as
% separate functions. In our case, we do not have nonlinear constraints. 
%--------------------------------------------------------------------------
function [c,ceq] = nonlcon(x)

    c = [];
    ceq = [];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program End