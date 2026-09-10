function [d,H_i] = PeriodFinder_n(P,i)

% input: stochastic matrix P of size(n,n), i state of P 
% output: d = period of d(i) 
%         H_i = set of states communicating with i 
n = size(P,1);
if (i > n | i < 1)
    error('invalid index i')
end

% finds powers of P up to n, P,P^2,...,P^n
R = PowersTransitionP(P);

% finds the partition of X into communication classes
C = DividingCommunicatingClasses(R); 

% finds number of communication classes
m = size(C,1);
d = 0; % initialize period to be 0 


% find communication class of i (H_i)
for j = 1:m
    if (C(j,i)==1)
    H_i = C(j,:); 
    break;
    end
end

n_i = sum(H_i); %number of elements in the class

% extract sub-matrix P_i sub transition matrix of H_i and all its powers 
indices = find(H_i);
R_i = cell(1,n_i);

for j = 1:n_i
    if (d==1)
        break
    end

R_i{j} = R{j}(indices,indices); % extract (P_i)^j

% update period
    if (diag(R_i{j}) ~= 0)
        d = gcd(d,j);
    end

end

end

function [C] = DividingCommunicatingClasses(R)
% input : R array of powers of P 
% output C = matrix whose rows are the communication classes
n = length(R);
C = zeros(n,n);
inserted = zeros(1,n); 
for i = 1:n

   
    if (inserted == ones(1,n))
        break;
    end
    if (inserted(i)==1)
        continue;
    end
    % 

    for j = 1:n
    if (i==j)
        C(i,i) = 1;
        inserted(i) = 1;
    end
    
    if ( i_communicates_j(R,i,j,n))
        C(i,j) = 1;
        inserted(j) = 1; 
    end
    end
end

C(all(C == 0, 2), :) = [];


end


function [flag] = i_communicates_j(R,i,j,n)
% this function checks if i comunicates with j 
% by checking if pji ^k >0 and pij^h >0 
% n size of P 
flag_ij = 0;
flag_ji = 0;
    for k = 1:n
        if (R{k}(i,j) > 0 & flag_ij ==0)
            flag_ij = 1; 
        end
        if (R{k}(j,i) > 0 & flag_ji ==0)
            flag_ji = 1; 
        end
        
    end
    flag = flag_ij*flag_ji;

end


function [R] = PowersTransitionP(P)
    % This function takes a stochastic matrix P 
    % and returns an array R such that:
    % R(1) = P, R(2) = P^2, ..., R(n) = P^n
if (IsStochastic(P))
    n = size(P,1);
    % Preallocate the cell array to hold the powers of P
    R = cell(1, n);
    
    % Initialize R(1) to be P
    R{1} = P;
    
    % Compute the powers of P iteratively
    for k = 2:n
        R{k} = R{k-1} * P;  % R(k) is P^(k) = P^(k-1) * P
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




