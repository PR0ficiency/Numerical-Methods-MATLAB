% ********************************
%  Divided Difference Calculator
% ********************************
% *
% *    Written by: Craig Medlin
% * Last Modified: March 9, 2016

function [ DD ] = divDiff(X, Y)
%divDiff This function calculates the divided differences and returns
% a table containing the results. 
%   INPUTS (required)
%         X - A vector of size n containing x values
%         Y - A vector of size n containing corresponding y values
%

% Input validation
if nargin < 2
    error('DivDiff requires 2 inputs!');
end % if statement

% Initialize Variables
n = size(X,1);
DD = zeros(n);
DD(:,1) = Y;

% Loop for each column
for j = 2 : n 
    
    % Loop for each row
    for i = 1:(n-j+1)
        
        % Store result
        DD(i,j) = (DD(i+1,j-1) - DD(i,j-1)) / (X(i+j-1) - X(i));
        
    end % i loop
    
end % j loop

end % function
