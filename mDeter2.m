function [ detA ] = mDeter2( A )
%M_DETER2 finds the determinant of an input square matrix.
% Improved version of mDeter
% Written by Craig Medlin - CC_BY 2016

rows = size(A,1);
cols = size(A,2);
detA = 0;
swaps = 0;

% Verify Square Matrix
if rows == cols
    n = rows;
    
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
    
    detA = (-1)^swaps;
    for i=1:n
        detA = detA*A(i,i);
    end
    
else
    fprintf('ERROR! mDeter() only accepts a square matrix argument.');
end % if

end % function mDeter

