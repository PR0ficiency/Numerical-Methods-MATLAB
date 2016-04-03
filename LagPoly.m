% *****************************************
%     Lagrangian Polynomial Interpolation
% *****************************************
% *
% *    Written by: Craig Medlin
% * Last Modified: March 10, 2016

function [ L ] = LagPoly(X, Y, xi)
%LagPoly This function interpolates a value for an unknown point
% given a set of x points and corresponding f(x) values using 
% the Lagrangian polynomial method.
%
%   INPUTS (required)
%         X - A set of n x values
%         Y - A set of n corresponding y values
%        xi - X value to interpolate for
% 
%   OUTPUTS 
%        L - Interpolated value at xi
%
%
%   EXAMPLE
%       X = [1.1 2.2 3.3 4.4];
%       Y = [6 12 3 6];
%       L = LabPoly(X,Y,3.1);
%

% Input validation
narginchk(3,inf);

% Initialization
L = 0; % Interpolated value
n = length(X);

% Loop for each point
for k = 1:n
    
    % Initialize product value
    P = 1;
    
    % Loop for each point
    for i = 1:n
        % Check if i is equal to k
        if i ~= k
            P = P * (xi - X(i))/(X(k)-X(i));
        end % If
    end % For
    
    % Update interpolation value
    L = L + P * Y(k);
    
end % For
        
