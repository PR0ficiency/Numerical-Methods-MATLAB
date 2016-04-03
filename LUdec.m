function [ L U ] = LUdec( A )
%LUDEC performs LU decomposition on an input matrix A.
% This function is a varient of the function GuassElim.m
% Written by Craig Medlin - CC_BY 2016

% Get size of input matrix
    n = size(A,1);
    swaps = 0;
    
    for j = 1:(n-1)
        % Initialization
        pivot_element = abs(A(j,j));
        pivot(j) = j;
        pivot_temp_row = j;

        % Find row to pivot
        for i=(j+1):n
            tmp_el = abs(A(i,j));
            if tmp_el > pivot_element
                pivot_element = tmp_el;
                pivot_temp_row = i;
            end % if
        end % for i

        % Switch rows if necessary
        if pivot(j) ~= pivot_temp_row
            A([j pivot_temp_row],:)=A([pivot_temp_row j],:);
            swaps = swaps + 1;
        end % if

        % Record element multipliers
        for i = (j+1):n
            A(i,j) = A(i,j)/A(j,j);
        end % for i

        % Transform to upper triangular matrix
        for i = (j+1):n
            for k = (j+1):n
                A(i,k) = A(i,k) - A(i,j)*A(j,k);
            end % for k
        end % for i

    end % for j
    
    L = diag(ones(1,n));
    U = zeros(n);
    
    for i = 1:n
        L(i+1:end,i)=A(i+1:end,i);
        U(1:i,i)=A(1:i,i);
    end
    
end

