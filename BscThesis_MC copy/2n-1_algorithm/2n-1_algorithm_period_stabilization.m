function [d] = PeriodFinder_2n_1(P,i)

% input: stochastic matrix P of size(n,n), i state of P 
% output: d = period of d(i) 
 

n = size(P,1);
M = eye(n);
d = 0; %initialize period

if (i > n)
    error('invalid index i')
end

if (~IsStochastic(P))
    error('\nMatrix P is not Stochastic\n');
end

for k = 1:(2*n-1)
   if (d == 1)
       break;
   end

    M = M*P;
    if (M(i,i) > 0)
        d = gcd(d,k);

    end


end
end



function [flag] = IsStochastic(P)
    % This function checks if the matrix P is square and stochastic
    % A matrix is stochastic if:
    % 1. It is square (number of rows equals number of columns)
    % 2. All elements are non-negative
    % 3. The sum of each row is 1
    
    % Check if the matrix is square
    [rows, cols] = size(P);
    if rows ~= cols
        flag = false;
        disp('Matrix is not square.');
        return; % Exit the function early if the matrix is not square
    end
    
    % Check if all elements are non-negative
    if any(P(:) < 0)
        flag = false;
        disp('Matrix contains negative elements.');
        return; % Exit the function early if negative elements are found
    end
    
    % Check if each row sums to 1
    rowSums = sum(P, 2); % Summing along the rows (second dimension)
    
    if all(abs(rowSums - 1) < 1e-10) % Allowing for small numerical tolerance
        flag = true;
    else
        flag = false;
        disp('Row sums are not equal to 1.');
    end
end






