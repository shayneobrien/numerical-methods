%% Algorithms 6.4, 6.5, 6.6

% Name: Shayne O'Brien
% Course: MATH 345 (Dr. Haddad)
% Due Date: Friday, 11/11/16 by 11:59 pm
% Content: Algorithm 6.4 used to solve 6.5 #7d  #8d
%          Algorithm 6.5 used to solve 6.6 #7c  #8c
%          Algorithm 6.6 used to solve 6.6 #7d  #8d

%% Preparation
clear all; % clear all 
format long; % set format to long

%% Algorithm 6.4: lU Factorization
% Modify Algorithm 6.4 (lU-decomp) to be able to solve a linear system. 
% Use it to solve: 6.5 #7d, 8d (see data file for this, comment out ones
% that are not desired. All results are included in the diary file in the
% zip folder.

% To factor the nxn matrix A = [a(i,j)] into the product of the
% lower-triangular matrix l = [l(i,j)] and the upper-triangular matrix U =
% [u(i,j)]; that is, A = lU, where the main diagonal of either l or U
% consists of all ones:

% INPUT: dimension n; the entries a(i,j), 1<=i, j<=n of A; the diagnology
% l(1,1) = ... = l(n,n) = 1 or the diagonoal u(1,1) = ... = u(n,n) = 1 of
% U.

% OUTPUT: the entires l(i,j), 1<=j<=i, 1<=i<=n of l and the entries u(i,j),
% i<=j<=n, 1 <=i<=n of U.

Section6AlgsData_OBrien; % load in the data
l = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]; %the diagonal of the lower triangular matrix all set to 1
u(1,1) = a(1,1)./l(1,1); % sumelect l(1,1) and u(1,1) satisfying l(1,1)*u(1,1) == a(1,1)  
                         % i.e., set the first value of u

sum = 0; % initialize sum. This will help us do summations later on

if l(1,1)*u(1,1) == 0 % If factorization is not possible:
    disp('Factorization impossible') % display error
else % Otherwise, we proceed with the rest of the algorithm:
    for j = 2:n % For entries 2 to n:
        u(1,j) = a(1,j)./l(1,1); %creates the first row of U 
        l(j,1) = a(j,1)./u(1,1); %creates the first column of l
    end % end for loop
    for i = 2:n-1 % For entries 2 to n-1:
        for k = 1:i-1 % For entries 1 to n-1:
            sum = sum + (l(i,k).*u(k,i)); % compute summation needed to calculate u(i,i)
        end % end for loop
        u(i,i) = (a(i,i) - sum)./l(i,i); % compute u(i,i) since we already know l(i,i)
        sum = 0; % set sum back to 0
        if l(i,i)*u(i,i) == 0 % If factorization is not possible:
            disp('Factorization impossible') % display error
        else % Otherwise, proceed:
            for j = i+1:n % For entries i+1 to n:
                for k = 1:i-1 % For entries 1 to i-1:
                    sum = sum + l(i,k).*u(k,j); % compute summation needed for u(i,j)
                end % end for loop
                u(i,j) = (1./l(i,i)).*(a(i,j) - sum); % calculate the ith column of l
                sum = 0; % set sum back to 0
                for k = 1:i-1 % For the first entry until the i-1 th entry
                    sum = sum + l(j,k).*u(k,i); % compute sum
                end % end for loop
                l(j,i) = (1./u(i,i)).*(a(j,i) - sum); % compute l(j,i)
                sum = 0; % reset sum to 0
            end % end for loop
        end % end if/else loop
    end
    sum = 0; % reset sum to 0;
    for k = 1:n-1 % For entries 1 to n-1:
            sum = sum + (l(n,k).*u(k,n)); % compute sum for finding u(n,n) term
    end % end for loop
    u(n,n) = (a(n,n) - sum)./l(n,n); % compute last value of diagnol of u
end % end loop

% We now have l and U, so we will display them
l
u

% Now we will solve the linear system
y = zeros(1,n); % initialize y
y(1)=b(1)./l(1,1);  % we will use forward solving first to solve ly=b. 
                    % the first y value is equal to the first b value since we have ones on the diagonal of our l
sum=0; % reinitialize sum as 0

for i = 2:n % For entries 2 to n:
    for j=1:i-1 % For entries 1 to i-1:  
                %calculate the summation term  that to help us find the rest of y
                % it will multiply the known y-values with their corresponding coefficients 
                % to use later to subtract from b to find the next y value
        sum = sum + l(i,j).*y(j); % compute sum  
    end % end for loop
    y(i) = (b(i) - sum)./l(i,i); % calculate the rest of the y-values
    sum = 0; % set sum back to 0
end % end for loop

% Now we back solve to solve Ux=y
x(n)=y(n)./u(n,n); % compute the n-th entry of x
sum = 0; % set sum back to 0
for i = n-1:-1:1; % Going backwards from the n-1 entry to 1:
    for j=i+1:n; % For entries i+1 to n:
        sum = sum + u(i,j).*x(j); %calculate the summation term that will help us in finding the other values of x
    end % end for loop
    x(i) = (y(i) - sum)./u(i,i); % calculate the rest of the x-values
    sum = 0; % set sum back to 0
end

%output x with the final answer
x

%% Algorithm 6.5: lDlt Factorization
% Modify Algorithm 6.5 (lDl^t) to be able to solve a linear system and 
% use it to solve 6.6 #7c, 8c (see data file for this, comment out ones
% that are not desired). All results are included in the diary file in the
% zip folder.

% To factor the positive definite nxn matrix A into the form lDlt, where l
% is a lower triangular matrix with 1s along the diagonal and D is a
% diagonal matrix with positive entries on the diagonal:

% INPUT: the dimension n; entries a(i,j) for 1<=i, j<=n of A
% OUTPUT: the entires l(i,j) for 1<=j<i and 1<=i<=n of l, and d(i) for
% 1<=i<=n of D.

Section6AlgsData_OBrien; % load in data
l = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]; %the diagonal of the lower triangular matrix all set to 1
d = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]; %d is a diagonal matrix with positive entries
sum = 0; %initialize sum at 0
d(1,1) = a(1,1); % set the first entry of d as equal to the first entry of a

% We will now compute the LDLt factorization
for i = 1:n % For entries 1 to n:
    for j = 1:i-1 % For entries 1 to i-1:
        v(j) = l(i,j).*d(j,j); % compute the j-th entry of v
        sum = sum + l(i,j).*v(j); % compute the sum to be used for calculating d(i,i)
    end % end for loop
    d(i,i) = a(i,i) - sum; % compute the diagonal of d
    sum = 0; % reset sum to 0
    for j = i+1:n % For entries i+1 to n:
        for k = 1:i-1 % For entries 1 to i-1:
            sum = sum + l(j,k).*v(k); % compute sum to be used for calculating l(j,i)
        end % end for loop
        l(j,i) = (a(j,i) - sum)./d(i,i); % compute l(j,i)
        sum=0; % reset sum back to 0
    end % end for loop
end % end for loop
sum=0; % reset sum to 0

% We now have l and d, so we will output them:
l % output l
d % output d

% Now we will solve the linear system by first solving ly=b
y(1) = b(1); % let the first entry of y be the first entry in b
for i = 2:n % For entries 2 to n:
    for j = 1:i-1 % For entries 1 to i-1:
        sum=sum + l(i,j).*y(j); % compute sum to be used to calculate y(i)
    end % end for loop
    y(i) = b(i) - sum; % compute y(i)
    sum=0; % reset sum
end % end for loop

% Now we will solve the system Dz=y
for i=1:n % For entries 1 to n:
    z(i) = y(i)./d(i,i); % compute z(i)
end % end for loop

% Finally, we will now solve Ltx = z
x(n) = z(n); % compute the n-th entry of x
for i = n-1:-1:1 % From entry n-1 backwards to 1:
    for j = i+1:n % For entries i+1 to n:
        sum = sum + l(j,i).*x(j); % compute sum to be used to calculate x(i)
    end % end for loop
    x(i) = z(i) - sum; % compute x(i)
    sum = 0; % reset sum to 0
end % end for loop

% Display final result x
x

%% Algorithm 6.6: Cholesky (llt)

% To factor the positive definite nxn matrix A into llt , where l is lower
% triangular.

% INPUT: the dimension n; entries a(i,j) for 1<=i, j<=n of A.
% OUTPUT: the entries l(i,j) for 1<=j<=i and 1<=i<=n of l. (The entries of
% U=lt are u(i,j) = l(j,i) for i<=j<=n and 1<=i<=n.

Section6AlgsData_OBrien; % load in the data
l = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]; %the diagonal of the lower triangular matrix all set to 1
sum = 0; % initialize sum to 0

l(1,1) = sqrt(a(1,1)); % compute l(1,1)
for j = 2:n % For entries 2 to n:
    l(j,1) = a(j,1)./l(1,1); % compute the columns of l
end % end for loop
for  i = 2:n-1 % For entries 2 to n-1:
    for k = 1:i-1 % For entries 1 to i-1: 
        sum = sum + (l(i,k))^2; % Compute sum to be used for l(i,i)
    end % end for loop
    l(i,i) = (a(i,i)-sum)^(1/2); % compute diagonal of l
    sum = 0; % set sum to 0
    for j = i+1:n % For entries i+1 to n:
        for k = 1:i-1 % For entries 1 to i-1:
            sum = sum + l(j,k).*l(i,k); % compute sum to calculate l(j,i)
        end % end for loop
        l(j,i) = (a(j,i)-sum)./l(i,i); % compute l(j,i)
        sum = 0; % reset sum to 0
    end % end for loop
end % end for loop
for k = 1:n-1 % For entries 1 to n-1:
    sum = sum + (l(n,k))^2; % compute sum to calculate diagonal of l
end % end for loop
l(n,n) = (a(n,n)-sum)^(1/2); % compute l(n,n)
sum = 0; % reset sum to 0
for i=1:n % For entries 1 to n:
    for j=i:n % For entries i to n:
        lT(i,j)=l(j,i); % for finding l-transpose
    end % end for loop
end % end for loop

% Output the l and l-transpose matrices
l 
lT

% Now we solve the linear system
y(1) = b(1)./l(1,1); % compute y(1)
for i = 2:n % For entries 2 to n:
    for j = 1:i-1; % For entries 1 to i-1: 
        sum = sum + l(i,j).*y(j); % compute sum for calculating y(i)
    end % end for loop
    y(i) = (b(i) - sum)./l(i,i); % compute y(i)
    sum = 0; % reset sum to 0
end % end for loop
x(n) = y(n)./l(n,n); % compute x(n)
for i = n-1:-1:1 % For entries n-1 backwards to 1:
    for j = i+1:n % For entries i+1 to n:
        sum = sum + l(j,i).*x(j); % compute sum used to calculate x(i)
    end % end for loop
    x(i) = (y(i) - sum)./l(i,i); % compute x(i)
    sum = 0; % reset sum to 0
end % end for loop

% Output the solution
x
