% this function will computes the LU factorization of A
% input: A; outputs: L, U

function [L, U] = gaussLU(A)
    % check if A is square
    [m, n] = size(A);
    if m ~= n
        error('Input matirx A should be square');
    end

    for k = 1 : n-1
        for i = k+1 : n
            piv = A(i,k) / A(k,k);
            if isnan(piv)
                error('zero pivot encountered');
            end
            A(i, k+1:n) = A(i, k+1:n) - piv * A(k, k+1:n); % r_curr = r_curr - piv * r_prev
            A(i, k) = piv; % restore pivs
        end % inner for loop end
    end % outer for loop end
    % get L and U
    L = tril(A, -1) + eye(n);
    U = triu(A);
end % end of gaussLU