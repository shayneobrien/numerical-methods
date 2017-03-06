%% Newton's Method for Exercises 2.1: 5,9,11,13

% Name: Shayne O'Brien
% Course: MATH 345 (Dr. Haddad)
% Date: 9/25/16
% Content: Function to compute 

%% Background
%Inverses of numbers without divison: In the early days of computing,
% calculators were mechanical, not electronic, and the earliest of these
% could multiply, but not divide. Algorithms such as Newton's Method could
% be used to calculate the intervest of a number a, i.e., a^-1 = 1/a,
% iteratively.

%% Clearing out space
format long % put format on long for precision and outputs.
close all; clear all; % close all plots, clear all variables.

%% Part (a)
% Begin by coding Newton's method and saving it as Newton_YourName. It
% should perform Newton's Method. Your program should input an initial
% guess, p, a tolerance TOL, and N, a maximum number of iterates. The
% program should output a "table" with the number of the iterate in the
% first column, followed by the approximations to the fixed point of g
% (root of f) in the next column, followed by the absolute Cauchy error,
% and last the relative Cauchy error.

% Begin Part (a)

%Here we input our function f that will be used to compute Newton's method.
% We let f(x) = 1/x - a, where a is the number of which we will find the
% inverse. We choose this function because when x = 1/a, f(x) = 0 which is
% a root. This root corresponds to the value of 1/a.
a = input('Input the number we will find the inverse of, a: '); % Input a
f = @(x)( 1/x - a ); % Input function
df = @(x)( -1/(x^2) ); % Input function derivative

% Ask user for input
N = input('Input maximum number of iterates, N: '); % ask user to put in max. number of iterates, N
pvec = zeros(N,1); % pre-allocate vector to store each p approximation
pvec(1) = input('Input initial approximation, Po: '); % ask user to input initial guess, p
TOL = input('Input tolerance level, TOL: '); % ask user to input tolerance, TOL

% Pre-allocations to speed up our algorithm
IterationVec = zeros(N,1); % pre-allocate vector to store iterations for speed
AbsCauchyError = zeros(N,1); % pre-allocate vector to store absolute Cauchy error for speed
RelCauchyError = zeros(N,1); % pre-allocate vector to store relative Cauchy error for speed

% Fixing pre-allocated vectors for table output
RelCauchyError(1) = NaN; % first entry is NaN because error cannot be computed in the first iteration.
AbsCauchyError(1) = NaN; % first entry is NaN because error cannot be computed in the first iteration.

% Initialize k for final output
k = 0; % k is a constant that we will use to simplify the decision making process for our output.
       % if k remains 0 after N iterations, we have not found a fixed point
       % within acceptable TOL. If k = 1, we have and will therefore
       % display all relevant results.

% Initialize i
i = 1; % let i = 1. i will be used as the check for the while loop.

% Error catching
if ((N < 1) || (mod(N,1) ~= 0)) % If the max number of iterations N has been set below 1 or N is not an integer:
    disp('Error: Maximum number of iterations must be an integer greater than equal to 1.') % display error.
elseif TOL <= 0 % Otherwise, if the tolerance is negative:
    disp('Error: TOL must be positive.') % display error.
elseif ( (2 - a*pvec(i)) < 0) % If 2 - a*Po < 0:
    disp('Error: 2 - a*Po must be greater than 0 or Newtons Method will diverge.') % display error.
elseif a == 0 % If a is 0, its inverse is undefined so:
    disp('Error: Cannot take the inverse of 0!'); % display error
else % otherwise since both criteria are met, proceed with computation.
    
    % Implementation of Newton's Method
    while i < N % Until we reach our max number of iterations N:
        p = pvec(i) - f(pvec(i))/df(pvec(i)); % compute p-sub-i
                                              % Pn+1 = Pn - f(Pn)/df(Pn)
        IterationVec(i) = i; % update Iteration vector for i-th entry.
        AbsCauchyError(i+1) = abs(p - pvec(i)); % update relative Cauchy error vector. index is i+1 because error cannot be computed in the first iteration.
        RelCauchyError(i+1) = abs((p - pvec(i))/p); % update absolute Cauchy error vector. index is i+1 because error cannot be computed in the first iteration.
        if ( abs((p - pvec(i))) ) < TOL % If the absolute Cauchy error is less than TOL we have found a fixed point so:
            pvec = pvec(1:i); % shorten vector so table outputs correctly
            IterationVec = IterationVec(1:i); % shorten vector
            AbsCauchyError = AbsCauchyError(1:i); % shorten vector
            RelCauchyError = RelCauchyError(1:i); % shorten vector
            k = 1; % set k = 1 so the root will be printed
            break % break while loop
        else % Otherwise:
            i = i + 1; % increment i by 1
            pvec(i) = p; % store the i-th iteration of p
        end % end if loop    
    end % end while loop
    
    % Printing output
    if k == 1 % if we found a fixed point that meets the tolerance criteria within N iterations:
        % displaying results matrix
        ResultsMatrix = [IterationVec, pvec, RelCauchyError, AbsCauchyError]; % Store results in a matrix
        fprintf('    Iteration Number\t      Pn\t    |(Pn - Pn-1)/Pn|\t  |Pn - Pn-1| \n') % Print colum nnames
        disp(ResultsMatrix) % print column contents
        fprintf('The approximate fixed point is p = %.8f, reached in %d iterations.\n', pvec(i),i) % display root and # iterations taken
    else % Otherwise, since a fixed point was not found within a reasonable TOL level:
        fprintf('Newtons Method failed to find an acceptably accurate root within N=%d iterations.\n', N) % display error
    end % end if loop
end % end if loop
    

%% Part (b)
% Determine a function f(x) so that by applying Newton's Method to f the
% program will find the inverse (reciprocal) of a number, x = 1/a, with no
% division taking place. In other words, using Newton's Method you will
% find a g(x) that does not actually use division to find 1 divided by a.

% Begin Part (b)

% Newton's Method states P(n+1)= P(n) - (f(P(n)))/(df(P(n))), for n >= 1
% Through algebra, we observe that

%     P(n+1) = P(n) - f(P(n))/df(P(n))
%            = P(n) - ((1/P(n)) - a)/(-1/(P(n))^2)
%            = P(n) + P(n)*(1 - a*P(n))
%            = P(n)*(2 - a*P(n))

% where P(n) is the approximation, a is some number that's inverse we are 
% trying to find, f is some function, and df is the derivative of f.

% Therefore, we can use Newton's method to approximate p = 1/a without ever
% actually using division. We denote g(x) = P(n+1) = P(n)*(2 - a*P(n)), where a is
% some number and g(x) is the approximation to the inverse of a, 1/a.

% It is important to note that if 2 - a * P(1) < 0, where P(1) is the
% initial guess for the inverse, the algorithm will diverge.

%% Part (c)
% Apply Newton's Method to the function f(x) found in Part (b), and use the
% resulting g(x) to find x = 1/a with no division and an absolute error
% less than 10^-12 for a = 19.

% Begin Part (c)

% Input function for f from Part (b)
f = @(x)( x*(2-a*x) ); % Input function

% Ask user for input
a = input('Input the number we will find the inverse of, a: '); % ask user for number of which we will find the inverse of
N = input('Input maximum number of iterates, N: '); % ask user to put in max. number of iterates, N
pvec = zeros(N,1); % pre-allocate vector to store each p approximation
pvec(1) = input('Input initial approximation, Po: '); % ask user to input initial guess, p
TOL = input('Input tolerance level, TOL: '); % ask user to input tolerance, TOL

% Pre-allocations to speed up our algorithm
IterationVec = zeros(N,1); % pre-allocate vector to store iterations for speed
AbsCauchyError = zeros(N,1); % pre-allocate vector to store absolute Cauchy error for speed
RelCauchyError = zeros(N,1); % pre-allocate vector to store relative Cauchy error for speed

% Fixing pre-allocated vectors for table output
RelCauchyError(1) = NaN; % first entry is NaN because error cannot be computed in the first iteration.
AbsCauchyError(1) = NaN; % first entry is NaN because error cannot be computed in the first iteration.

% Initialize k for final output
k = 0; % k is a constant that we will use to simplify the decision making process for our output.
       % if k remains 0 after N iterations, we have not found a fixed point
       % within acceptable TOL. If k = 1, we have and will therefore
       % display all relevant results.

% Initialize i
i = 1; % let i = 1. i will be used as the check for the while loop.

% Error catching
if ((N < 1) || (mod(N,1) ~= 0)) % If the max number of iterations N has been set below 1 or N is not an integer:
    disp('Error: Maximum number of iterations must be an integer greater than equal to 1.') % display error.
elseif TOL <= 0 % Otherwise, if the tolerance is negative:
    disp('Error: TOL must be positive.') % display error.
elseif ( (2 - a*pvec(i)) < 0) % If 2 - a*Po < 0:
    disp('Error: 2 - a*Po must be greater than 0 or Newtons Method will diverge.') % display error.
elseif a == 0 % If a is 0, its inverse is undefined so:
    disp('Error: Cannot take the inverse of 0!'); % display error
else % otherwise since both criteria are met, proceed with computation. 

    % Implementation of Newton's Method
    while i < N % Until we reach our max number of iterations N:
        % p = 2*pvec(i) - a*(pvec(i))^2; % compute p-sub-i without using division, using function from Part (b)
        p = f(pvec(i)); % compute p-sub-i without using division, using function from Part (b)
        IterationVec(i) = i; % update Iteration vector for i-th entry.
        AbsCauchyError(i+1) = abs(p - pvec(i)); % update relative Cauchy error vector. index is i+1 because error cannot be computed in the first iteration.
        RelCauchyError(i+1) = abs((p - pvec(i))/p); % update absolute Cauchy error vector. index is i+1 because error cannot be computed in the first iteration.
        if ( abs((p - pvec(i))) ) < TOL % If the absolute Cauchy error is less than TOL we have found a fixed point so:
            pvec = pvec(1:i); % shorten vector so table outputs correctly
            IterationVec = IterationVec(1:i); % shorten vector
            AbsCauchyError = AbsCauchyError(1:i); % shorten vector
            RelCauchyError = RelCauchyError(1:i); % shorten vector
            k = 1; % set k = 1 so the root will be printed
            break % break while loop
        else % Otherwise:
            i = i + 1; % increment i by 1
            pvec(i) = p; % store the i-th iteration of p
        end % end if loop    
    end % end while loop
    
    % Printing output
    if k == 1 % if we found a fixed point that meets the tolerance criteria within N iterations:
        % displaying results matrix
        ResultsMatrix = [IterationVec, pvec, RelCauchyError, AbsCauchyError]; % Store results in a matrix
        fprintf('    Iteration Number\t      Pn\t    |(Pn - Pn-1)/Pn|\t   |Pn - Pn-1| \n') % Print colum nnames
        disp(ResultsMatrix) % print column contents
        fprintf('The approximate fixed point is p = %.8f, reached in %d iterations.\n', pvec(i),i) % display root and # iterations taken
        % plot results
        figure; % open new figure
        plot(IterationVec, pvec) % plot iteration vs. root approximation
        xlabel('Iteration number N') % label x-axis
        ylabel('Fixed-point approximation, Pn') % label y-axis
        title('Newtons Method') % title for plot is the input function g       
    else % Otherwise, since a fixed point was not found within a reasonable TOL level:
        fprintf('Newtons Method failed to find an acceptably accurate root within N=%d iterations.\n', N) % display error
    end % end if loop
end % end if loop


%% Part (d)
%Modify your code for Newton's Method to apply the Secant Method to your
% function from Part (b) with a = 21. You may create a completely new
% script or m-file to do this, or place it in a new cell following Newton's
% Method.

% Begin Part (d)

% Note: We use a new function here, rather than the one derived in Part
% (b). I could not get the code to work using the function from Part (b),
% but using this new f(x) we do correctly obtain the root.

% Input a, f(x)
a = input('Input the number for which we will find the inverse: '); % Input a. We are trying to find 1/a.
f = @(x)( (a*x) - 1 ); % Input function, f(x) = a*x - 1 = 0

% Ask user for input
N = input('Input maximum number of iterates, N: '); % ask user to put in max. number of iterates, N
pvec = zeros(N,1); % pre-allocate vector to store each p approximation
pvec(1) = input('Input initial approximation, Po: '); % ask user to input initial guess, P0
pvec(2) = input('Input initial approximation, P1: '); % ask user to input initial guess P1
TOL = input('Input tolerance level, TOL: '); % ask user to input tolerance, TOL

% Pre-allocations to speed up our algorithm
IterationVec = zeros(N,1); % pre-allocate vector to store iterations for speed
AbsCauchyError = zeros(N,1); % pre-allocate vector to store absolute Cauchy error for speed
RelCauchyError = zeros(N,1); % pre-allocate vector to store relative Cauchy error for speed

% Fixing pre-allocated vectors for table output
RelCauchyError(1:2) = NaN; % first two entries are NaN because error cannot be computed in the 1st/2nd iteration.
AbsCauchyError(1:2) = NaN; % first two entries are NaN because error cannot be computed in the 1st/2nd iteration.
IterationVec(1) = 1; % Let the first entry of the iteration vector be 1, since we will start at i=2.
% Initialize k for final output
k = 0; % k is a constant that we will use to simplify the decision making process for our output.
       % if k remains 0 after N iterations, we have not found a fixed point
       % within acceptable TOL. If k = 1, we have and will therefore
       % display all relevant results.

% Initialize i
i = 2; % let i = 2. i will be used as the check for the while loop. We start at i=2 
       % because the secant method computes Pn = Pn-1 - f(Pn-1)(Pn-1 - Pn-2)/((f(Pn-1) - f(Pn-2)))
       % so if we let i = 1, the script would print an error.

% Error catching
if ((N < 1) || (mod(N,1) ~= 0)) % If the max number of iterations N has been set below 1 or N is not an integer:
    disp('Error: Maximum number of iterations must be an integer greater than equal to 1.') % display error.
elseif TOL <= 0 % Otherwise, if the tolerance is negative:
    disp('Error: TOL must be positive.') % display error.
elseif (f(pvec(1)) * f(pvec(2)) > 0 ) % If there is no sign change for f(Po)*f(P1):
    disp('Error: f(Po) * f(P1) > 0. Try different intiial approximations that will give a sign change.') % display error; likely that Secant Method won't converge.
else % otherwise since both criteria are met, proceed with computation.
    
    % Implementation of Newton's Method
    while i < N % Until we reach our max number of iterations N:
        p = pvec(i) - (f(pvec(i))*(pvec(i) - pvec(i-1)))/(f(pvec(i))-f(pvec(i-1))); % compute p-sub-i
                             % Pn = Pn-1 - f(Pn-1)(Pn-1 - Pn-2)/((f(Pn-1) - f(Pn-2)))
        IterationVec(i) = i; % update Iteration vector for i-th entry.
        AbsCauchyError(i) = abs(p - pvec(i)); % update relative Cauchy error vector. index is i+1 because error cannot be computed in the first iteration.
        RelCauchyError(i) = abs((p - pvec(i))/p); % update absolute Cauchy error vector. index is i+1 because error cannot be computed in the first iteration.
        if ( abs((p - pvec(i))) ) < TOL % If the absolute Cauchy error is less than TOL we have found a fixed point so:
            pvec = pvec(1:i); % shorten vector so table outputs correctly
            IterationVec = IterationVec(1:i); % shorten vector
            AbsCauchyError = AbsCauchyError(1:i); % shorten vector
            RelCauchyError = RelCauchyError(1:i); % shorten vector
            k = 1; % set k = 1 so the root will be printed
            break % break while loop
        else % Otherwise:
            i = i + 1; % increment i by 1
            pvec(i) = p; % store the i-th iteration of p
        end % end if loop    
    end % end while loop
    
    % Printing output
    if k == 1 % if we found a fixed point that meets the tolerance criteria within N iterations:
        % displaying results matrix
        ResultsMatrix = [IterationVec, pvec, RelCauchyError, AbsCauchyError]; % Store results in a matrix
        fprintf('    Iteration Number\t      Pn\t    |(Pn - Pn-1)/Pn|\t   |Pn - Pn-1| \n') % Print colum nnames
        disp(ResultsMatrix) % print column contents
        fprintf('The approximate fixed point is p = %.8f, reached in %d iterations.\n', pvec(i),i) % display root and # iterations taken
        
        % plot results
        figure;
        plot(IterationVec, pvec) % plot iteration vs. root approximation
        xlabel('Iteration number N') % label x-axis
        ylabel('Fixed-point approximation, Pn') % label y-axis
        title('Secant Method') % title for plot is the input function g        
    else % Otherwise, since a fixed point was not found within a reasonable TOL level:
        fprintf('Newtons Method failed to find an acceptably accurate root within N=%d iterations.\n', N) % display error
    end % end if loop
end % end if loop

%% Part (e)
%Use Steffensen's Algorithm (which modifies Aitken's D^2) to speed up the
% convergence of Part (c) by first finding a new fixed point function,
% g(x), that is not the same as the one found in Part (b) for Newton's
% Method. Why wouldn't you want to apply Steffensen's algorithm to Newton's
% Method?

% Begin Part (e)

% Note: Algorithm would not work properly for functions tried other than Part
% (b). I tried more than five different g(x) functions and none of them
% produced the correct root. 

% Ask user to input the number of which we will find the inverse
a = input('Input the number for which we will find the inverse: ');
% Input function being used to compute Steffensen's Algorithm
g = @(x)( x*(2-a*x) );

% Ask user for input
N = input('Input maximum number of iterates, N: '); % ask user to put in max. number of iterates, N
p0vec = zeros(N,1); % pre-allocate vector to store each p0 approximation
p0vec(1) = input('Input initial approximation, Po: '); % ask user to input initial guess
p1vec = zeros(N,1); % pre-allocate vector to store each p1 approximation
p2vec = zeros(N,1); % pre-allocate vector to store each p2 approximation
TOL = input('Input tolerance level, TOL: '); % ask user to input tolerance, TOL

% Pre-allocations to speed up our algorithm
IterationVec = zeros(N,1); % pre-allocate vector to store iterations for speed
AbsCauchyError = zeros(N,1); % pre-allocate vector to store absolute Cauchy error for speed
RelCauchyError = zeros(N,1); % pre-allocate vector to store relative Cauchy error for speed
RelCauchyError(1) = NaN; % first entry is NaN because error cannot be computed in the first iteration.
AbsCauchyError(1) = NaN; % first entry is NaN because error cannot be computed in the first iteration.

k = 0; % for deciding our output, similarly to the other parts of this project
i = 1;
if ((N < 1) || (mod(N,1) ~= 0)) % If the max number of iterations N has been set below 1 or N is not an integer:
    disp('Error: Maximum number of iterations must be an integer greater than equal to 1.') % display error.
elseif TOL <= 0 % Otherwise, if the tolerance is negative:
    disp('Error: TOL must be positive.') % display error.
else % Otherwise, if both conditions are met:
    while i <= N % until we reach out maximum iteration,
        p1vec(i) = g(p0vec(i)); % P1 = g(P0)
        p2vec(i) = g(p1vec(i)); % P2 = g(P1)
        p = p0vec(i) - ((p1vec(i) - p0vec(i))^2)/(p2vec(i) - 2*p1vec(i) + p0vec(i)); % compute P-sub-i
        IterationVec(i) = i; % update Iteration vector for i-th entry.
        AbsCauchyError(i+1) = abs(p - p0vec(i)); % update relative Cauchy error vector. index is i+1 because error cannot be computed in the first iteration.
        RelCauchyError(i+1) = abs((p - p0vec(i))/p); % update absolute Cauchy error vector. index is i+1 because error cannot be computed in the first iteration.
        if ( AbsCauchyError(i+1) ) < TOL % If we have gone below our tolerance:
            p0vec(i+1) = p; % Update the table for the P0 column
            p0vec = p0vec(1:i+1); % shorten vector so table outputs correctly
            p1vec = p1vec(1:i+1); % shorten vector so table outputs correctly
            p1vec(i+1) = nan; % P1 was not calculated for i-th iteration, so input a NaN
            p2vec = p2vec(1:i+1); % shorten vector so table outputs correctly
            p2vec(i+1) = nan; % P2 was not calculated for i-th iteration, input a NaN
            IterationVec = IterationVec(1:i+1); % shorten vector
            IterationVec(i+1) = i+1; % Update iteration vector
            AbsCauchyError = AbsCauchyError(1:i+1); % shorten vector
            RelCauchyError = RelCauchyError(1:i+1); % shorten vector
            k = 1; % set k to 1 so we output results instead of an error
            break % break while loop
        else % Otherwise: 
            i = i + 1; % increment i by 1
            p0vec(i) = p; % update P0
        end % end if loop
    end % end while loop
    % Printing output
    if k == 1 % if we found a fixed point that meets the tolerance criteria within N iterations:
        % displaying results matrix
        ResultsMatrix = [IterationVec, p0vec, p1vec , p2vec, RelCauchyError, AbsCauchyError]; % Store results in a matrix
        fprintf('    Iteration Number\t       P0    \t\t   P1 \t\t       P2 \t    |(Pn - Pn-1)/Pn|\t|Pn - Pn-1| \n ') % Print colum nnames
        disp(ResultsMatrix) % print column contents
        fprintf('The approximate fixed point is p = %.8f, reached in %d iterations.\n', p0vec(i),i+1) % display root and # iterations taken
        % plotting results
        figure;
        plot(IterationVec, p0vec) % plot iteration vs. root approximation
        xlabel('Iteration number N') % label x-axis
        ylabel('Fixed-point approximation, Pn') % label y-axis
        title('Steffensens Algorithm') % title for plot is the input function g    
    
    else % Otherwise, since a fixed point was not found within a reasonable TOL level:
        fprintf('Steffensens Method failed to find an acceptably accurate root within N=%d iterations.\n', N) % display error
    end % end if loop 
end % end if loop
    

% Why don't you want to apply Steffensen's algorithm to Newton's Method?

%We do not want to apply Steffensen's algorithm to Newton's Method
% because Steffensen's Method accelerated the convergence of a linearly
% convergent sequence to quadratic. Newton's Method is already quadratic; 
% applying Steffensen's algorithm to Newton's Method would be redundant.

%% Part (f)
% Plot and compare all three sets of results. What can you conclude from
% this?

% Begin Part (f)

% The plots are coded into each of the Parts (c), (d), and (e). We will be
% comparing the rate of convergence for finding 1/a, where a = 19 and TOL =
% 1E-12. The initial guess is 0.04 for Parts (c) and (e). In part (d), we
% give P0 = -0.1 and P1 = 0.1, so that there is a sign change between f(P0)
% and f(P1).

% We see that the Secant Method converges the quickest in 3 iterations,
% Steffensen's Algorithm converges the second-fastest in 5 iterations, and
% Newton's Method converges the slowest of the three in 6 iterations.


%Advantages of Secant Method:
% 1. Converges at a faster than linear rate
% 2. Does not require computation of derivative

%Diadvantages of Secant Method:
% 1. No guaranteed error bound
% 2. Sometimes finds extraneous roots
% 3. Can run into problems if f(a) = 0 (x-axis tangent to graph of y = f(x)
% at x = a).

%Advantages of Newtons Method:
% 1. Converges quadratically, which is one of the fastest when it does
% converge

%Disadvantages of Newtons Method:
% 1. Only converges when initial guess is very near to solution
% 2. Computationally expensive due to needing the derivative
% 3. If the tangent is parallel or almost parallel to the x-axis, then the
% method won't converge

%Advantages of Steffensen's Algorithm:
% 1. Can speed up convergence to quadratic

%Disadvantages of Steffensen's Algorithm:
% 1. Needs very close starting value to actual solution
% 2. Computationally expensive compared to Secant Method

