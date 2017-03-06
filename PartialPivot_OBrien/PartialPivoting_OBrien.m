%% Algorithm 6.2 and 6.3: Gaussian Elimination with Pivoting

% Name: Shayne O'Brien
% Course: MATH 345 (Dr. Haddad)
% Due Date: Wednesday, 11/2/16 by 11:30 am
% Content: Algorithm 6.2,6.3 code, computations for all problems

%% a. 
% Please hand in from 6.1 #1d, 10 by hand.  
% Please hand in from 6.2 #9b,10c,11b  by hand.  

% These were done by hand and turned in on Wednesday, 11/2/16.

%% b. 
%Please code Algorithm 6.2 and use your code to solve 6.2 #13b, 15b, 14c. 
% You will have to make modifications to perform rounding vs. chopping, etc.  
% (Chopping means you cut off the remainder (2-digit chopping means 1.29 becomes 1.2.  
%  Rounding means it becomes 1.3.) 
%Be sure to test your code on at least 2 examples from the book or class to make sure that it works.  
% You should calculate the absolute error and the relative error in each case.

% Each of these three problems is separated into a different cell.

%% #13b: Use Gaussian elimination (with partial pivoting) to solve the 
% following linear systems, and compare the approximations to the actual 
% solution:

% 3.03*x1 - 12.1*x2 + 14*x3 = -119,
% -3.03*x1 + 12.1*x2 - 7*x3 = 120,
% 6.11*x1 - 14.2*x2 + 21*x3 = -139
% Actual solution [0,10,1/7].

clear all; % clear all variables

actual = [0;10;1/7]; % actual solution
A = [3.03, -12.1, 14; -3.03, 12.1, -7; 6.11, -14.2, 21]; % system of equations
B = [-119; 120; -139]; % what they are equal to

% We will code Algorithm 6.2: Gaussian Elimination with Partial Pivoting,
% which will solve the n x n linear system

% E1: a(11)x(1) + a(12)x(2) + ... + a(1n)x(n) = a(1,n+1)
% E2: a(21)x(1) + a(22)x(2) + ... + a(2n)x(n) = a(2,n+1)
%                              .
%                              .
%                              .
% En: a(n1)x(1) + a(n2)x(2) + ... + a(nn)x(n) = a(n,n+1)

% INPUT: number of unknowns and equations n; augmented matrix A = [a(ij)],
% where 1 <= i <= n and 1 <= j <= n+1

% OUTPUT: solution x1,...,xn or message that linear system has no unique
% solution.

matrix = [A, B]; % augment matrix A and column B to form a matrix representation of our linear system
n = rank(A); % n is the rank of our augmented matrix A
NROW = zeros(n,1); % let NROW be a column of zeros

for i = 1:n % For all entries
    NROW(i) = i; % initialize row pointer
end % end for loop

for i = 1:n-1 % For the first entry until the penultimate entry:
    maxval = 0; % initialize maxval as a variable equal to 0
    p = 0; % initialize p as a variable equal to 0
    for j = i:n % For every entry between current entry i and the last entry n, we will find the max val in the column:
        if abs(matrix(NROW(j),i)) > maxval % If the current maxval is less than the value being checked
            maxval = abs(matrix(NROW(j),i)); % update maxval
            p = j; % update the p
        end % end if loop
    end % end for loop
    if matrix(NROW(p),i) == 0 % If this entry is zero, we have an equation equal to zero so
        disp('No unique solution exists.') % print no unique solution exists
        break % break the loop
    else % Otherwise, we row interchange:
        if NROW(i) ~= NROW(p) % If they aren't already equal:
            NCOPY = NROW(i); % copy NROW(i) to a variable
            NROW(i) = NROW(p); % let entry i be entry p
            NROW(p) = NCOPY; % let entry p be entry i
        end % end if loop
        for k = (i+1):n % For every entry between i+1 and the last entry n, we will do Gaussian elimination:
            m(NROW(k),i) = matrix(NROW(k),i)/matrix(NROW(i),i); %calculate the coefficient needed to row reduce
            for l = i:n+1 % for all entries from i to n+1
                matrix(NROW(k),l) = matrix(NROW(k),l) - m(NROW(k),i)*matrix(NROW(i),l); % row reduce
            end % end for loop
        end % end for loop
    end % end if loop
end % end for loop

x = zeros(n,1); % let x be a row of 0s of length n
x(n) = matrix(NROW(n),n+1)/matrix(NROW(n),n); % compute last x(n)

for i = n-1:-1:1 % for the rest of the x's
    sum = 0; % let sum be 0 
    for j = i+1:n % for all entries from i+1 to n
        sum = sum + matrix(NROW(i),j)*x(j); % compute sum
    end % end for loop
    x(i) = (matrix(NROW(i),n+1) - sum)/matrix(NROW(i),i); % compute x(i)
end % end while loop

abserr = norm(actual-x); % compute Absolute Error
relerr = norm(actual-x)/norm(x); % compute Relative Error

% display results
fprintf('Approximation: [%.3f; %.3f; %.3f]\nAbsolute error: %f \nRelative error: %f\n',x,abserr,relerr)

% Printout should be
% Approximation: [0.000,10.000,0.143]
% Absolute error: 0.000000
% Relative error: 0.000000

% These errors are very low; this is a very good approximation.

%% 6.2 #15b
% #15b: Repeat Exercise 9b using Gaussian elimination with partial pivoting
% and three-digit rounding arithmetic.

% The system of equations is given in the previous cell. We will now edit the code
% from the that cell to account for three-digit rounding arithmetic.

clear all; % clear all variables

actual = [0;10;1/7]; % this is the actual answer
A = [3.03, -12.1, 14; -3.03, 12.1, -7; 6.11, -14.2, 21]; % system of equations
B = [-119; 120; -139]; % what they are equal to

matrix = [A, B]; % augment matrix A and column B to form a matrix representation of our linear system
n = rank(A); % n is the rank of our augmented matrix A
NROW = zeros(n,1); % let NROW be a column of zeros

for i = 1:n % For all entries
    NROW(i) = i; % initialize row pointer
end % end for loop

for i = 1:n-1 % For the first entry until the penultimate entry:
    maxval = 0; % initialize maxval as a variable equal to 0
    p = 0; % initialize p as a variable equal to 0
    for j = i:n % For every entry between current entry i and the last entry n, we will find the max val in the column:
        if abs(matrix(NROW(j),i)) > maxval % If the current maxval is less than the value being checked
            maxval = abs(matrix(NROW(j),i)); % update maxval We do not have to round here because maxval is never used for computation.
            p = j; % update the p
        end % end if loop
    end % end for loop
    if matrix(NROW(p),i) == 0 % If this entry is zero, we have an equation equal to zero so
        disp('No unique solution exists.') % print no unique solution exists
        break % break the loop
    else % Otherwise, we row interchange:
        if NROW(i) ~= NROW(p) % If they aren't already equal:
            NCOPY = NROW(i); % copy NROW(i) to a variable
            NROW(i) = NROW(p); % let entry i be entry p
            NROW(p) = NCOPY; % let entry p be entry i
        end % end if loop
        for k = (i+1):n % For every entry between i+1 and the last entry n, we will do Gaussian elimination:
            m(NROW(k),i) = round(matrix(NROW(k),i),3,'significant')/round(matrix(NROW(i),i),3,'significant'); %calculate the coefficient needed to row reduce
            m(NROW(k),i) = round(m(NROW(k),i),3,'significant'); % round to 3 sig figs
            for l = i:n+1 % for all entries from i to n+1
                matrix(NROW(k),l) = round(matrix(NROW(k),l),3,'significant') - round(m(NROW(k),i),3,'significant')*round(matrix(NROW(i),l),3,'significant'); % row reduce
                matrix(NROW(k),l) = round(matrix(NROW(k),l),3,'significant'); % round to 3 sig figs
            end % end for loop
        end % end for loop
    end % end if loop
end % end for loop

x = zeros(n,1); % let x be a row of 0s of length n
x(n) = round(matrix(NROW(n),n+1),3,'significant')/round(matrix(NROW(n),n),3,'significant'); % compute last x(n)
x(n) = round(x(n),3,'significant'); % round to 3 sig figs

for i = n-1:-1:1 % for the rest of the x's
    sum = 0; % let sum be 0 
    for j = i+1:n % for all entries from i+1 to n
        sum = round(sum,3,'significant') + round(matrix(NROW(i),j)*x(j),3,'significant'); % compute sum
        sum = round(sum,3,'significant'); % round
    end % end for loop
    x(i) = (round(matrix(NROW(i),n+1),3,'significant') - round(sum,3,'significant'))/round(matrix(NROW(i),i),3,'significant'); % compute x(i)
    x(i) = round(x(i),3,'significant'); % round
end % end while loop

abserr = norm(actual-x); % compute Absolute Error
relerr = norm(actual-x)/norm(x); % compute Relative Error

% display results
fprintf('\nApproximation: [%.3f; %.3f; %.3f].\n Absolute Error: %f\n Relative error: %f\n',x,abserr,relerr)

% Print out should be
% Approximation: [0.000; 10.000; 0.143].
% Absolute Error: 0.000143
% relative error: 0.000014

% We can see that by rounding using three-digit arithmetic, the error goes
% up. Still, the approximation for this problem is acceptable since the error
% is not very high.

%% 6.2 #14c
% #14c: Use Gaussian elimination with partial pivoting and three-digit
% chopping arithmetic to solve the following linear systems and compare the
% approximations to the actual solution.

% 1.19*x1 + 2.11*x2 - 100*x3 + x4 = 1.12
% 14.2*x1 - 0.122*x2 + 12.2*x3 - x4 = 3.44
%           100*x2 - 99.9*x3 + x4 = 2.15
% 15.3*x1 + 0.110*x2 - 13.1*x3 - x4 = 4.16

% Actual solution is [0.176, 0.0126, -0.0206, -1.18].

clear all; % clear all variables

actual = [0.176; 0.0126; -0.0206; -1.18];
A = [1.19, 2.11, -100, 1; 14.2, -0.122, 12.2, -1; 0, 100, -99.9, 1; 15.3, 0.110, -13.1, -1];
B = [1.12; 3.44; 2.15; 4.16];

matrix = [A, B]; % augment matrix A and column B to form a matrix representation of our linear system
n = rank(A); % n is the rank of our augmented matrix A
NROW = zeros(n,1); % let NROW be a column of zeros

for i = 1:n % For all entries
    NROW(i) = i; % initialize row pointer
end % end for loop

for i = 1:n-1 % For the first entry until the penultimate entry:
    maxval = 0; % initialize maxval as a variable equal to 0
    p = 0; % initialize p as a variable equal to 0
    for j = i:n % For every entry between current entry i and the last entry n, we will find the max val in the column:
        if abs(matrix(NROW(j),i)) > maxval % If the current maxval is less than the value being checked
            maxval = abs(matrix(NROW(j),i)); % update maxval
            p = j; % update the p
        end % end if loop
    end % end for loop
    if matrix(NROW(p),i) == 0 % If this entry is zero, we have an equation equal to zero so
        disp('No unique solution exists.') % print no unique solution exists
        break % break the loop
    else % Otherwise, we row interchange:
        if NROW(i) ~= NROW(p) % If they aren't already equal:
            NCOPY = NROW(i); % copy NROW(i) to a variable
            NROW(i) = NROW(p); % let entry i be entry p
            NROW(p) = NCOPY; % let entry p be entry i
        end % end if loop
        for k = (i+1):n % For every entry between i+1 and the last entry n, we will do Gaussian elimination:
            m(NROW(k),i) = chop(matrix(NROW(k),i),3)/chop(matrix(NROW(i),i),3); %calculate the coefficient needed to row reduce
            m(NROW(k),i) = chop(m(NROW(k),i),3); % chop to 3 sig figs
            for l = i:n+1 % for all entries from i to n+1
                matrix(NROW(k),l) = chop(matrix(NROW(k),l),3) - chop(m(NROW(k),i),3)*chop(matrix(NROW(i),l),3); % row reduce
                matrix(NROW(k),l) = chop(matrix(NROW(k),l),3); % chop to 3 sig figs
            end % end for loop
        end % end for loop
    end % end if loop
end % end for loop

x = zeros(n,1); % let x be a row of 0s of length n
x(n) = chop(matrix(NROW(n),n+1),3)/chop(matrix(NROW(n),n),3); % compute last x(n)
x(n) = chop(x(n),3); % chop to 3 sig figs

for i = n-1:-1:1 % for the rest of the x's
    sum = 0; % let sum be 0 
    for j = i+1:n % for all entries from i+1 to n
        sum = chop(sum,3) + chop(matrix(NROW(i),j)*x(j),3); % compute sum
        sum = chop(sum,3); % chop
    end % end for loop
    x(i) = (chop(matrix(NROW(i),n+1),3) - chop(sum,3))/chop(matrix(NROW(i),i),3); % compute x(i)
    x(i) = chop(x(i),3); % chop
end % end while loop

abserr = norm(actual-x); % compute Absolute Error
relerr = norm(actual-x)/norm(x); % compute Relative Error

% display results
fprintf('\nApproximation: [%.3f; %.3f; %.3f; %.3f].\n Absolute Error: %f\n Relative error: %f\n',x,abserr,relerr)

% Print out should be
% Approximation: [0.178; 0.013; -0.020; -1.160].
% Absolute Error: 0.020101
% Relative error: 0.017124

% The error for this approximation is ~2% absolute, ~1.7% relative. This is
% decent, but not great.


%% Additional checks for correctness

% I ran the above code also on the following two problems:
% 1. "Illustration" from page 376:
% A = [30, 591400; 5.291, -6.130];
% B = [591700; 46.78];
% actual = [10, 1]

% I got:

% partial pivoting (regular)
% x = [10, 1]
% abserr = 0;
% relerr = 0;

% partial pivoting (3-digit rounding)
% x = [10, 1]
% abserr = 0;
% relerr = 0;

% partial pivoting (3-digit chopping)
% x = [10, 1]
% abserr = 0;
% relerr = 0;

% This is a very good approximation.

% 2. Example 10d that begins on page 376:
% A = [pi, sqrt(2), -1, 1; exp(1), -1, 1, 2; 1, 1, -sqrt(3), 1; -1, -1, 1, -sqrt(5)]
% B = [0;1;2;3];
% actual = [1.35; -4.68; -4.03; -1.66].

% I got:

% partial pivoting (regular)
% Approximation: [1.349; -4.678; -4.033; -1.656638]
% Absolute error: 0.004902 
% Relative error: 0.001

% partial pivoting (3-digit rounding)
% Approximation: [1.360; -4.720; -4.070; -1.650].
% Absolute Error: 0.058310
% Relative error: 0.008849

% partial pivoting (3-digit chopping)
% Approximation: [1.360; -4.720; -4.070; -1.650].
% Absolute Error: 0.058310
% Relative error: 0.008849

% Chopping and rounding reduce the accuracy of the approximation. Regularly,
% the approximation is good. The absolute error for chopping and rounding
% is a bit high.


%% c. 
%Please modify your code to perform Algorithm 6.3 and use it to solve 
% 6.2 #17b, 19b, 18c.  
% Be sure to test your code on at least 2 examples from the book or class 
% to make sure that it works.  You should calculate the absolute error and 
% the relative error in each case.

% The only steps that differ from those of Algorithm 6.2 are the first,
% second, and third step.

%% 6.2 #17b
% Repeat Exercise 9b using Gaussian elimination with scaled partial pivoting

% The system of equations is given in the previous cell. I have edited the
% code to do scaled partial pivoting.

clear all; % clear all variables

actual = [0;10;1/7]; % this is the actual answer
A = [3.03, -12.1, 14; -3.03, 12.1, -7; 6.11, -14.2, 21]; % system of equations
B = [-119; 120; -139]; % what they are equal to

matrix = [A, B]; % augment matrix A and column B to form a matrix representation of our linear system
n = rank(A); % n is the rank of our augmented matrix A
NROW = zeros(n,1); % let NROW be a column of zeros
s = zeros(n,1); % intialize vector to store s
maxval = 0; % initialize maxval as a variable equal to 0
p = 0; % initialize p as a variable equal to 0


s = max(abs(A'))'; % this makes a vector of the maximums of each row. this
                   % is a more efficient way to find the maximums than to
                   % use a nested for loop
for i = 1:n % For all entries
    NROW(i) = i; % initialize row pointer
    if s(i) == 0 % if there is a zero in any of the maximums
        disp('No unique solution exists.') % no unique solution exists
        break % break the loop
    end % end for loop
end % end for loop

for i = 1:n-1 % For the first entry until the penultimate entry:
    maxval = 0; % reset maxval as a variable equal to 0
    p = 0; % reset p as a variable equal to 0
    for j = i:n % For every entry between current entry i and the last entry n, we will find the max val in the column:
        if abs(matrix(NROW(j),i))/s(NROW(j)) > maxval % If the current maxval is less than the value being checked
            maxval = abs(matrix(NROW(j),i))/s(NROW(j)); % update maxval
            p = j; % update the p
        end % end if loop
    end % end for loop
    if matrix(NROW(p),i) == 0 % If this entry is zero, we have an equation equal to zero so
        disp('No unique solution exists.') % print no unique solution exists
        break % break the loop
    else % Otherwise, we row interchange:
        if NROW(i) ~= NROW(p) % If they aren't already equal:
            NCOPY = NROW(i); % copy NROW(i) to a variable
            NROW(i) = NROW(p); % let entry i be entry p
            NROW(p) = NCOPY; % let entry p be entry i
        end % end if loop
        for k = (i+1):n % For every entry between i+1 and the last entry n, we will do Gaussian elimination:
            m(NROW(k),i) = matrix(NROW(k),i)/matrix(NROW(i),i); %calculate the coefficient needed to row reduce
            for l = i:n+1 % for all entries from i to n+1
                matrix(NROW(k),l) = matrix(NROW(k),l) - m(NROW(k),i)*matrix(NROW(i),l); % row reduce
            end % end for loop
        end % end for loop
    end % end if loop
end % end for loop

x = zeros(n,1); % let x be a row of 0s of length n
x(n) = matrix(NROW(n),n+1)/matrix(NROW(n),n); % compute last x(n)

for i = n-1:-1:1 % for the rest of the x's
    sum = 0; % let sum be 0 
    for j = i+1:n % for all entries from i+1 to n
        sum = sum + matrix(NROW(i),j)*x(j); % compute sum
    end % end for loop
    x(i) = (matrix(NROW(i),n+1) - sum)/matrix(NROW(i),i); % compute x(i)
end % end while loop

abserr = norm(actual-x); % compute Absolute Error
relerr = norm(actual-x)/norm(x); % compute Relative Error

% display results
fprintf('Approximation: [%.3f; %.3f; %.3f]\nAbsolute error: %f \nRelative error: %f\n',x,abserr,relerr)

% Print out should be
% Approximation: [0.000; 10.000; 0.143].
% Absolute Error: 0.000000
% relative error: 0.000000

% This is a very good approximation.


%% 6.2 #19b 
% Repeat Exercise 9b using Gaussian elimination with scaled partial pivoting
% and three-digit rounding arithmetic

% The system of equations is given in the previous cell. We will now edit the code
% from the that cell to account for three-digit rounding arithmetic with scaled partial pivoting.

clear all; % clear all variables

actual = [0;10;1/7]; % this is the actual answer
A = [3.03, -12.1, 14; -3.03, 12.1, -7; 6.11, -14.2, 21]; % system of equations
B = [-119; 120; -139]; % what they are equal to

matrix = [A, B]; % augment matrix A and column B to form a matrix representation of our linear system
n = rank(A); % n is the rank of our augmented matrix A
NROW = zeros(n,1); % let NROW be a column of zeros

s = round(max(abs(A'))',3,'significant'); % this makes a vector of the maximums of each row. this
                   % is a more efficient way to find the maximums than to
                   % use a nested for loop
for i = 1:n % For all entries
    NROW(i) = i; % initialize row pointer
    if s(i) == 0 % if there is a zero in any of the maximums
        disp('No unique solution exists.') % no unique solution exists
        break % break the loop
    end % end for loop
end % end for loop

for i = 1:n-1 % For the first entry until the penultimate entry:
    maxval = 0; % initialize maxval as a variable equal to 0
    p = 0; % initialize p as a variable equal to 0
    for j = i:n % For every entry between current entry i and the last entry n, we will find the max val in the column:
        if abs(matrix(NROW(j),i))/s(NROW(j)) > maxval % If the current maxval is less than the value being checked
            maxval = abs(matrix(NROW(j),i))/s(NROW(j)); % update maxval. We do not have to round here because maxval is never used for computation.
            p = j; % update the p
        end % end if loop
    end % end for loop
    if matrix(NROW(p),i) == 0 % If this entry is zero, we have an equation equal to zero so
        disp('No unique solution exists.') % print no unique solution exists
        break % break the loop
    else % Otherwise, we row interchange:
        if NROW(i) ~= NROW(p) % If they aren't already equal:
            NCOPY = NROW(i); % copy NROW(i) to a variable
            NROW(i) = NROW(p); % let entry i be entry p
            NROW(p) = NCOPY; % let entry p be entry i
        end % end if loop
        for k = (i+1):n % For every entry between i+1 and the last entry n, we will do Gaussian elimination:
            m(NROW(k),i) = round(matrix(NROW(k),i),3,'significant')/round(matrix(NROW(i),i),3,'significant'); %calculate the coefficient needed to row reduce
            m(NROW(k),i) = round(m(NROW(k),i),3,'significant'); % round to 3 sig figs
            for l = i:n+1 % for all entries from i to n+1
                matrix(NROW(k),l) = round(matrix(NROW(k),l),3,'significant') - round(m(NROW(k),i),3,'significant')*round(matrix(NROW(i),l),3,'significant'); % row reduce
                matrix(NROW(k),l) = round(matrix(NROW(k),l),3,'significant'); % round to 3 sig figs
            end % end for loop
        end % end for loop
    end % end if loop
end % end for loop

x = zeros(n,1); % let x be a row of 0s of length n
x(n) = round(matrix(NROW(n),n+1),3,'significant')/round(matrix(NROW(n),n),3,'significant'); % compute last x(n)
x(n) = round(x(n),3,'significant'); % round to 3 sig figs

for i = n-1:-1:1 % for the rest of the x's
    sum = 0; % let sum be 0 
    for j = i+1:n % for all entries from i+1 to n
        sum = round(sum,3,'significant') + round(matrix(NROW(i),j)*x(j),3,'significant'); % compute sum
        sum = round(sum,3,'significant'); % round
    end % end for loop
    x(i) = (round(matrix(NROW(i),n+1),3,'significant') - round(sum,3,'significant'))/round(matrix(NROW(i),i),3,'significant'); % compute x(i)
    x(i) = round(x(i),3,'significant'); % round
end % end while loop

abserr = norm(actual-x); % compute Absolute Error
relerr = norm(actual-x)/norm(x); % compute Relative Error

% display results
fprintf('\nApproximation: [%.3f; %.3f; %.3f]\n Absolute Error: %f\n Relative error: %f\n',x,abserr,relerr)

% Print out should be
% Approximation: [0.000; 10.000; 0.143].
% Absolute Error: 0.000143
% relative error: 0.000014

% We can see that by rounding using three-digit arithmetic, the error goes
% up. Still, the approximation for this problem is acceptable since the error
% is not very high. The error is the same for scaled partial pivoting as it
% is for partial pivoting.

%% 6.2 #18c. 
% Repeat Exercise 10c using Gaussian elimination with scaled partial
% pivoting and chopping.

% The system of equations is given in 6.2 #14c.

% 1.19*x1 + 2.11*x2 - 100*x3 + x4 = 1.12
% 14.2*x1 - 0.122*x2 + 12.2*x3 - x4 = 3.44
%           100*x2 - 99.9*x3 + x4 = 2.15
% 15.3*x1 + 0.110*x2 - 13.1*x3 - x4 = 4.16

% Actual solution is [0.176, 0.0126, -0.0206, -1.18].

clear all; % clear all variables

actual = [0.176; 0.0126; -0.0206; -1.18];
A = [1.19, 2.11, -100, 1; 14.2, -0.122, 12.2, -1; 0, 100, -99.9, 1; 15.3, 0.110, -13.1, -1];
B = [1.12; 3.44; 2.15; 4.16];

matrix = [A, B]; % augment matrix A and column B to form a matrix representation of our linear system
n = rank(A); % n is the rank of our augmented matrix A
NROW = zeros(n,1); % let NROW be a column of zeros

s = chop(max(abs(A'))',3); % this makes a vector of the maximums of each row. this
                   % is a more efficient way to find the maximums than to
                   % use a nested for loop
for i = 1:n % For all entries
    NROW(i) = i; % initialize row pointer
    if s(i) == 0 % if there is a zero in any of the maximums
        disp('No unique solution exists.') % no unique solution exists
        break % break the loop
    end % end for loop
end

for i = 1:n-1 % For the first entry until the penultimate entry:
    maxval = 0; % initialize maxval as a variable equal to 0
    p = 0; % initialize p as a variable equal to 0
    for j = i:n % For every entry between current entry i and the last entry n, we will find the max val in the column:
        if abs(matrix(NROW(j),i))/s(NROW(j)) > maxval % If the current maxval is less than the value being checked
            maxval = abs(matrix(NROW(j),i))/s(NROW(j)); % update maxval. We do not have to round here because maxval is never used for computation.
            p = j; % update the p
        end % end if loop
    end % end for loop
    if matrix(NROW(p),i) == 0 % If this entry is zero, we have an equation equal to zero so
        disp('No unique solution exists.') % print no unique solution exists
        break % break the loop
    else % Otherwise, we row interchange:
        if NROW(i) ~= NROW(p) % If they aren't already equal:
            NCOPY = NROW(i); % copy NROW(i) to a variable
            NROW(i) = NROW(p); % let entry i be entry p
            NROW(p) = NCOPY; % let entry p be entry i
        end % end if loop
        for k = (i+1):n % For every entry between i+1 and the last entry n, we will do Gaussian elimination:
            m(NROW(k),i) = chop(matrix(NROW(k),i),3)/chop(matrix(NROW(i),i),3); %calculate the coefficient needed to row reduce
            m(NROW(k),i) = chop(m(NROW(k),i),3); % chop to 3 sig figs
            for l = i:n+1 % for all entries from i to n+1
                matrix(NROW(k),l) = chop(matrix(NROW(k),l),3) - chop(m(NROW(k),i),3)*chop(matrix(NROW(i),l),3); % row reduce
                matrix(NROW(k),l) = chop(matrix(NROW(k),l),3); % chop to 3 sig figs
            end % end for loop
        end % end for loop
    end % end if loop
end % end for loop

x = zeros(n,1); % let x be a row of 0s of length n
x(n) = chop(matrix(NROW(n),n+1),3)/chop(matrix(NROW(n),n),3); % compute last x(n)
x(n) = chop(x(n),3); % chop to 3 sig figs

for i = n-1:-1:1 % for the rest of the x's
    sum = 0; % let sum be 0 
    for j = i+1:n % for all entries from i+1 to n
        sum = chop(sum,3) + chop(matrix(NROW(i),j)*x(j),3); % compute sum
        sum = chop(sum,3); % chop
    end % end for loop
    x(i) = (chop(matrix(NROW(i),n+1),3) - chop(sum,3))/chop(matrix(NROW(i),i),3); % compute x(i)
    x(i) = chop(x(i),3); % chop
end % end while loop

abserr = norm(actual-x); % compute Absolute Error
relerr = norm(actual-x)/norm(x); % compute Relative Error

% display results
fprintf('\nApproximation: [%.3f; %.3f; %.3f; %.3f]\n Absolute Error: %f\n Relative error: %f\n',x,abserr,relerr)

% Print out should be
% Approximation: [0.179; 0.013; -0.020; -1.150].
% Absolute Error: 0.030151
% Relative error: 0.025901

% The error for this approximation is ~3% absolute, ~2.5% relative. This is
% decent, but not great. We would prefer error at least a magnitude of 10
% lower.


%% Additional checks for correctness

% I ran the above code also on the following two problems:

% 1. 10a on page 380:
% A = [58.9, 0.03; -6.10, 5.31];
% B = [59.2; 47];
% actual = [1; 10]

% I got:
% partial scaled pivoting (regular)
% Approximation: [1; 10]
% Absolute Error: 0.000000
% Relative error: 0.000000

% partial scaled pivoting (3-digit rounding)
% Approximation: [1; 10]
% Absolute Error: 0.000000
% Relative error: 0.000000

% partial scaled pivoting (3-digit chopping)
% Approximation: [1; 10]
% Absolute Error: 0.000000
% Relative error: 0.000000

% These are all very good approximations.

% 2. Example 10d that begins on page 376:
% A = [pi, sqrt(2), -1, 1; exp(1), -1, 1, 2; 1, 1, -sqrt(3), 1; -1, -1, 1, -sqrt(5)]
% B = [0;1;2;3];
% actual = [1.35; -4.68; -4.03; -1.66];

% I got:
% partial scaled pivoting (regular)
% Approximation: [1.349; -4.678; -4.033 -1.657]
% Absolute error: 0.004902 
% Relative error: 0.000750

% partial scaled pivoting (3-digit rounding)
% Approximation: [1.360; -4.720; -4.070; -1.650]
% Absolute Error: 0.058310
% Relative error: 0.008849

% partial scaled pivoting (3-digit chopping)
% Approximation: [1.360; -4.720; -4.070; -1.650]
% Absolute Error: 0.058310
% Relative error: 0.008849

% Chopping and rounding increases the error. The absolute error for chopping and
% rounding is a bit high, but regularly it is a decent estimation (~0.5%
% error).