function [ x, A2, b2, detA ] = GaussElim( A, b )
%GAUSS_ELIM performs Guassian Elimination to solve a matrix
%   Transforms a supplied matrix (A) to an upper triangular matrix to
%   simplify solving Ax = b. Pivots rows and performs back substitution and
%   scaling to solve for x. Returns modified  A and b matricies, as well as
%   solution x.
% Written by Craig Medlin - CC_BY 2016


    % Get size of input matrix
    rows = size(A,1);
    
    swaps = 0;
    detA = 0;
    
    for j = 1:(rows-1)
        % Initialization
        pivot_element = abs(A(j,j));
        pivot(j) = j;
        pivot_temp_row = j;

        % Find row to pivot
        for i=(j+1):rows
            tmp_el = abs(A(i,j));
            if tmp_el > pivot_element
                pivot_element = tmp_el;
                pivot_temp_row = i;
            end % if
        end % for i

        % Switch rows if necessary
        if pivot(j) ~= pivot_temp_row
            A([j pivot_temp_row],:)=A([pivot_temp_row j],:);
            b([j pivot_temp_row],:)=b([pivot_temp_row j],:);
            swaps = swaps + 1;
        end % if

        % Record element multipliers
        for i = (j+1):rows
            A(i,j) = A(i,j)/A(j,j);
        end % for i

        % Transform to upper triangular matrix
        for i = (j+1):rows
            for k = (j+1):rows
                A(i,k) = A(i,k) - A(i,j)*A(j,k);
            end % for k
            b(i) = b(i) - A(i,j)*b(j);
        end % for i

    end % for j

    % Perform back substitution
    x(rows) = b(rows)/A(rows,rows);
    for j = (rows-1):-1:1
        x(j) = b(j);
        
        for k = rows:-1:(j+1)
            x(j) = x(j) - x(k)*A(j,k);
        end % for k
        x(j) = x(j)/A(j,j);
        
    end % for j
    
    % Calculate Determinate
    detA = (-1)^swaps;
    for i = 1:rows
        detA = detA*A(i,i);
    end
    
    A2 = A;
    b2 = b;
    
end % function GaussElim

