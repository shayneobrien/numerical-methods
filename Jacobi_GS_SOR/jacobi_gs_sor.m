%% Algorithms 7.1 (Jacobi), 7.2 (Gauss-Seidel), and 7.3 (SOR)

% Name: Shayne O'Brien
% Course: MATH 345 (Dr. Haddad)
% Due Date: Friday, 11/28/16 by 11:59 pm
% Content: Jacobi, Gauss-Seidel, SOR algorithms

%% Instructions
% In addition to the usual things you should include in your programs:
% a. I would like for you to use the book?s algorithms for the Jacobi, Gauss-Seidel, and SOR.  
% b. I would also like for you to use the matrix versions (not in the book, but easy to program), and compare them.  
% c. You should also be checking results by comparing x(k) to x actual using relative error, or if not known, use the 
%    upper bounds given in the book. You should also look at the residual.  

% Solve:
% 7.3 #11c, d  (these all involve the same matrix)
% 7.3 #16
% 7.3 #18 
% 7.4 #10 (goes with 7.3 #18)
% 7.5 #5c

% The data for these are in the other file. Uncomment desired data and run
% code.

format long; % change format to long for greater precision

%% Algorithm 7.1 Jacobi Iterative

% COMPARED TO OTHER METHODS:
% Jacobi Iterative is an iterative technique to solve the nxn linear system Ax=b that starts 
% with an initial approximation vector x(0) to the solution x and generates
% a sequence of vectors {x(k)} (from k = 0:infinity) that converges to x. A
% residual vector is defined as r = b - Ax', where x' is an
% approximation to the solution of the linear system defined by A=bx. This
% is what we use to check the TOL.

% To solve Ax=b given an initial approximation x(0):

% INPUT: the number of equations and unknowns n; the entries a(i,j), 1<=i,
% j<=n of the matrix A; the entries b(i), 1<=i<=n of bl the entries XO(i),
% 1<=i<=n of XO = x(0); tolerance TOL; max number of iterations N.

% OUTPUT: the approximate solution x(1), ..., x(n) or a message that the
% number of iterations was exceeded.


clear all; % clear any existing variables
data; % uncomment the data you would like to use in the data file
TOL = input ('Enter to the COMMAND WINDOW, TOL,  the tolerance \n'); % TOL allows us to specify when we should stop with the iterations
N = input ('Enter to the COMMAND WINDOW, N,  the maximum number of iterations \n'); % set maximum number of iterations, so if it doesnt converge and never gets to limit the program will stop

k = 1; % let k be 1 for the first iteration
s = 0; % let s be 0. We will use this variable for summation functions
while k <= N % While k is less than or equal to our max iterations N:
    for i = 1:n % For all entries 1 to n: 
        for j = 1:n % For all entries 1 to n:
            if j ~= i % If j is not equal to i, since we exclude this term in our approximations:
                s = s + a(i,j).*XO(j); % compute s
            end % end if loop
        end % end for loop
        x(i) = (1./a(i,i)).*(-s + b(i)); % calculate the approximate ith iteration solution for vector x
        s = 0; % reset s back to zero
    end % end for loop
    if max( x' - XO ) < TOL % If our change between iterations is less than our tolerance:
        x % print x
        break; % break the while loop
    end % end if loop
    k = k+1; % Otherwise, our tolerance has not been reached, so increment iterative count by one to go through the loop again
    for i = 1:n % For all entries 1 to n:
        XO(i) = x(i)'; % set our x vector equal to XO for the next iteration
    end % end for loop
end % end while loop

fprintf('Maximum number of iterations exceeded \n') % output message to user
XO % print the final approximation we get either within our TOL or our max iterations

% Now we need to calculate the error
A = [a b]; % cocatenate the system of equations
R = rref(A);  % first we need the actual solution to the system, so we calculate the reduced row echelon form of A
for i = 1:n % For entries 1 to n:
    XA(i)=R(i,n+1); % compute the actual solution
end % end for loop

disp('Actual solution: ')
XA = XA' % transpose XA and print the actual solution
% and calculate relative error
ER = max(abs(XA-XO))./max(XA); % calculate error
fprintf('ER: %.6f (reached after %d iterations)\n', ER, k) % print error

%% Algorithm 7.2 Gauss-Seidel Iterative

% COMPARED TO OTHER METHODS:
%A possible improvement in Algorithm 7.1 can be seen by reconsidering our
% equations. The components of x(k-1) are used to compute the componenets
% of x(k). But for i > 1, the compoenent of x1(k),...,xi-1(k) of x(k) have
% already been computed. Using these instead of x(k-1), we can get a better
% approximation.

% To solve Ax=b given an initial approximation x(0):

% INPUT: the number of equations and unknowns n; the entries a(i,j), 1<=i,
% j<=n of the matrix A; the entries b(i), 1<=i<=n of bl the entries XO(i),
% 1<=i<=n of XO = x(0); tolerance TOL; max number of iterations N.

% OUTPUT: the approximate solution x(1), ..., x(n) or a message that the
% number of iterations was exceeded.

clear all; % clear any existing variables
Data_JGS_OBrien; % uncomment the data you would like to use in the data file
a = [1 0 -2; -.5 1 -.25; 1 -.5 1];
n = rank(a); %dimensions of our matrix
b = [.2; -1.425; 2];
XO = [0; 0; 0];
TOL = input ('Enter to the COMMAND WINDOW, TOL,  the tolerance \n'); % TOL allows us to specify when we should stop with the iterations
N = input ('Enter to the COMMAND WINDOW, N,  the maximum number of iterations \n'); % set maximum number of iterations, so if it doesnt converge and never gets to limit the program will stop

k = 1; % let k be 1 for the first iteration
s = 0; % let s be 0. We will use this variable for summation functions
p = 0; % let p be 0. We will use this variable for summation functions
while k <= N % While k is less than or equal to our max iterations N
    for i = 1:n % For entries 1 to n:
        for j = 1:i-1 % For entries 1 to i-1:
            s = s + a(i,j).*x(j); % compute the first summation
        end % end for loop
        for j = i+1:n % For entries i+1 to n:
            p = p + a(i,j).*XO(j); % compute second summation
        end % end for loop
        x(i) = (1./a(i,i)).*(-s - p + b(i)); % calculate the values in our approximated x vector
        s = 0; % reset s back to zero
        p = 0; % reset p back to zero
    end % end for loop
    if max(x' - XO) < TOL % If our change between iterations is less than our tolerance:
        x % print the approximation
        break; % break while loop
    end % end if loop
    k = k+1; % Otherwise, our tolerance has not been reached, so increment iterative count by one to go through the loop again
    for i = 1:n % For entries 1 to n:
        XO(i) = x(i)'; % set our x vector equal to XO for the next iteration
    end % end for loop
end % end while loop
fprintf('Maximum number of iterations exceeded \n') % output message to user
XO % print final approximation we get with our tolerance and max iterations

% Now we need to calculate the error
A = [a b]; % cocatenate the system of equations
R = rref(A);  % first we need the actual solution to the system, so we calculate the reduced row echelon form of A
for i = 1:n % For entries 1 to n:
    XA(i)=R(i,n+1); % compute the actual solution
end % end for loop

disp('Actual solution: ')
XA = XA' % transpose XA and print the actual solution
% and calculate relative error
ER = max(abs(XA-XO))./max(XA); % calculate error
fprintf('ER: %.6f (reached after %d iterations)\n', ER, k) % print error

%% Algorithm 7.3 SOR

% COMPARED TO OTHER METHODS:
% We are interested in choices of w with 1 < w, and these are called
% over-relaxation methods. They are used to accelerate the convergence for
% systems that are convergent by the Gauss-Seidel technique (if we let w =
% 1, SOR becomes equivalent to Gauss-Seidel). The methods are abbrieviated
% SOR, for Successive Over-Relaxation, and are particularly useful for
% solving the linear systems that occur in the numerical solution of
% certain partial-differential equations.

% To solve Ax=b given an initial approximation x(0):

% INPUT: the number of equations and unknowns n; the entries a(i,j), 1<=i,
% j<=n of the matrix A; the entries b(i), 1<=i<=n of bl the entries XO(i),
% 1<=i<=n of XO = x(0); the parameter w; tolerance TOL; max number of 
% iterations N.

% OUTPUT: the approximate solution x(1), ..., x(n) or a message that the
% number of iterations was exceeded.

clear all; % clear any existing variables
Data_JGS_OBrien; % uncomment the data you would like to use in the data file
TOL = input ('Enter to the COMMAND WINDOW, TOL,  the tolerance \n'); % TOL allows us to specify when we should stop with the iterations
N = input ('Enter to the COMMAND WINDOW, N,  the maximum number of iterations \n'); % set maximum number of iterations, so if it doesnt converge and never gets to limit the program will stop
w = input ('Enter to the COMMAND WINDOW, w,  the parameter \n'); % ask for our w parameter value 
% (remember, w < 1 is under-relaxation, w = 1 is Gauss-Seidel, 1 < w is over-relaxation.)

k = 1; % let k be 1 for the first iteration
s = 0; % let s be 0. We will use this variable for summation functions
p = 0; % let p be 0. We will use this variable for summation functions
while k <= N % While k is less than or equal to our max iterations N
    for i = 1:n % For entries 1 to n:
        for j = 1:i-1 % For entries 1 to i-1:
            s = s + a(i,j).*x(j); % compute first summation
        end % end for loop
        for j = i+1:n % For entries i+1 to n:
            p = p + a(i,j).*XO(j); % compute second summation
        end % end for loop
        x(i) = (1-w).*XO(i) + (1./a(i,i)).*(w.*(-s - p + b(i))); % calculate the values in our approximated x vector
        s = 0; % reset s back to zero
        p = 0; % reset p back to zero
    end % end for loop
    if max(x' - XO) < TOL % If our change between iterations is less than our tolerance:
        x % print the approximation
        break; % break while loop
    end % end if loop
    k = k+1; %if tolerance is not reached, we increase our iterative count by one and go through the loop again
    for i = 1:n % For entries 1 to n:
        XO(i) = x(i); % set our x vector equal to XO for the next iteration
    end % end for loop
end % end while loop
fprintf('Maximum number of iterations exceeded \n') % output message to user
XO % print final approximation we get with our tolerance and max iterations

% Now we need to calculate the error
A = [a b]; % cocatenate the system of equations
R = rref(A);  % first we need the actual solution to the system, so we calculate the reduced row echelon form of A
for i = 1:n % For entries 1 to n:
    XA(i)=R(i,n+1); % compute the actual solution
end % end for loop

disp('Actual solution: ')
XA = XA' % transpose XA and print the actual solution
% and calculate relative error
ER = max(abs(XA-XO))./max(XA); % calculate error
fprintf('ER: %.6f (reached after %d iterations)\n', ER, k) % print error
