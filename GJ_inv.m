function [ b ] = GJ_inv( A )
%GJ_INV Finds the inverse using the Gauss Jordan Method
% Written by Craig Medlin - CC_BY 2016

    % Get size of input matrix
    rows = size(A,1);
    b = eye(rows);
    
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
                b(i,k) = b(i,k) - A(i,j)*b(j,k);
            end % for k            
            
        end % for i
        
        for i = 1:(j-1)
            for k = (j+1):rows
                A(i,k) = A(i,k) - A(i,j)*A(j,k);
                b(i,k) = b(i,k) - A(i,j)*b(j,k);
            end % for k            
        end % for j

    end % for j
    
end

