% ********************************
%     Multi-variate Newton Solver
% ********************************
% *
% *    Written by: Craig Medlin
% * Last Modified: March 8, 2016

function [ x0, iterations, fnorm ] = mvNewton(A, b, eps, x0)
%mvNewton This function solves linear equations in the form Ax=b using 
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