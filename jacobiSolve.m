% ********************************
%     Jacobi Iterative Solver
% ********************************
% *
% *    Written by: Craig Medlin
% * Last Modified: March 8, 2016

function [ x0, iterations, fnorm ] = jacobiSolve(A, b, eps, x0)
%JACOBIITER This function solves linear equations in the form Ax=b using 
% the Jacobi Iterative mnethod until a certain error tolerance is met.
%
%   INPUTS (required)
%         A - A square system of n equations
%         b - A n x 1 system
%
%   INPUTS (optional)
%       eps - Error tolerance
%        x0 - Initial Guess of size n x 1
% 

% Initializations
n = size(A,1);
D = zeros(n);
O = zeros(n);
iterations = 0;

%% Input validation
if n ~= size(A,2)
    error('ERROR in jacobiSolve(): A matrix must be square!');
end  % IF Statement

% Check for minimum required arguments
if nargin < 2
    error('ERROR in jacobiSolve(): Not enough input arguments.');
end % IF Statement

% Check for error tolerance
if nargin < 3
    % Set default error tolerance
    eps = 10^-4;
end % IF Statement


%% Matrix Decomposition

% Loop once for each row
for i = 1:n
    % Extract diagonal elements
    D(i,i) = A(i,i);
    
    % Loop once for each column
    for j = 1:n
        if i~=j
            O(i,j) = A(i,j);
        end % IF Statement
    end % For each column
 
end % For each row
 

%% Iterate solution

% Calculate true solution
x_true = A\b;

% Set initial guess, if necessary
if nargin < 4
    x0 = D\b;
end

% Loop until tolerance met
while (norm(x_true - x0)) > eps
    
    % Calculate x_k
    x0 = D\(b - O*x0);
    
    % Increment iterations
    iterations = iterations + 1;
    
end % While loop

% Return final norm
fnorm = norm(x_true - x0);

end

