function [ detA ] = mDeter( A )
%M_DETER finds the determinant of an input square matrix.
%
% Written by Craig Medlin - CC_BY 2016

rows = size(A,1);
cols = size(A,2);
detA = 0;

% Verify Square Matrix
if rows == cols
    n = rows;
    
    if n > 2
        for i=1:n
            % Set addition or subtraction
            if mod(i,2) ~= 0
                j = 1;
            else
                j = -1;
            end % if
            
            % Columns selector 
            cols = linspace(1,n,n);
            cols = cols(cols~=i);
            
            % Create sub matrix
            tmp = A(2:end, cols);
            
            % Calculate determinant
            detA = detA + j * A(1,i)*mDeter(tmp);
        end % for
        
    else
        detA = A(1,1)*A(2,2) - A(1,2)*A(2,1);
    end % if
    
else
    fprintf('ERROR! mDeter() only accepts a square matrix argument.');
end % if

end % function mDeter

