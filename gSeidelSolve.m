% ********************************
%     Gauss-Seidel Iterative Solver
% ********************************
% *
% *    Written by: Craig Medlin
% * Last Modified: March 8, 2016

function [ x0, iterations, fnorm ] = gSeidelSolve(A, b, eps, x0)
%gSeidelSolve This function solves linear equations in the form Ax=b using 
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
iterations = 0;

%% Input validation
if n ~= size(A,2)
    error('ERROR in gSeidelSolve(): A matrix must be square!');
end  % IF Statement

% Check for minimum required arguments
if nargin < 2
    error('ERROR in gSeidelSolve(): Not enough input arguments.');
end % IF Statement

% Check for error tolerance
if nargin < 3
    % Set default error tolerance
    eps = 10^-4;
end % IF Statement

%% Iterate Solution

% Calculate true solution
x_true = A\b;

% Set initial guess, if necessary
if nargin < 4  
    D = zeros(n);
    % Loop once for each row
    for i = 1:n
        % Extract diagonal elements
        D(i,i) = A(i,i);
    end % For each row
 
    x0 = D\b;
end % IF Statement

% Loop until tolerance exceeded
while ( norm(x_true - x0) > eps )
    % For each row
    for i=1:n
        % Clear temp variable
        temp = 0;
        
        % For each column
        for j = 1:n
            
            % Make sure indicies do not match
            if j ~= i
                temp = temp + A(i,j)*x0(j);
            end % IF Statement
            
        end % For each column
        
        x0(i) = (b(i) - temp) / A(i,i);

    end % for each row
    
    % Increment iterations
    iterations = iterations + 1;
end % While

fnorm = norm(x_true - x0);

end % gSeidelSolve()

