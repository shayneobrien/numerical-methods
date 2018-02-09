%% Bisection Method for Exercises 2.1: 5,9,11,13

% Name: Shayne O'Brien
% Course: MATH 345 (Dr. Haddad)
% Date: 9/11/16
% Content: Function to compute the Bisection Method


%% This is the m-file that will perform the bisection method for various functions.
% Your code should prompt for the endpoints of the interval and tolerance.X
% It should verify a root exists and then prompt for new ones, if not.X
% You should output the results in tabular form, plot the sequences
% containing the endpoints and iterations (on the same plot).  
% You should be able to generate roots for many different functions.


%% Function to compute the Bisection Method
function bisect( f, a, b, TOL, N )

%Given the function f, left endpoint a, right endpoint b, tolerance threshold TOL, 
% and maximum number of iterations to be computed N, this % function will
% use the Bisection Method to find or approximate the root of function @f. 

%Usage note: it is necessary a function handle (@) before the input
% argument f. For example, to call the function x-2^(-x), the input for f
% would be @(x)x-2^(-x). Alternatively, we can create another .m function
% file and call the file name of that .m file using a function handle.

%Begin function contents:
format long % put format on long for precision and outputs.
k = 0; % this variable will be used to determine if we should output an error
       % or an answer later on in the function.
i = 1; % Initialize the first iteration by assigning i to the value 1.
FA = f(a); % Compute f(a), where f is some specified function.
FB = f(b); % Compute f(b). This will only be used to check that the f(a) * f(b) < 0.

% Each of the following will be used for the tabular output later on.
% By preallocating these vectors, the computation time is reduced greatly.
LeftEndPoint = zeros(N,1); % Preallocate vector for Left End Points (a)
RightEndPoint = zeros(N,1); % Preallocate vector for Right End Points (b)
Approximation = zeros(N,1); % Preallocate vector for Approximation (p)
FunctionValue = zeros(N,1); % Preallocate vector for Function Value (f(p))
ErrorBound = zeros(N,1); % Preallocate vector for Error Bound (f(p))

if b < a % If endpoint b is less than endpoint a,
    disp('Error: Endpoint a must be less than endpoint b.') % display error.
elseif ((N < 1) || (mod(N,1) ~= 0)) % If the max number of iterations N has been set below 1 or N is not an integer,
    disp('Error: Maximum number of iterations must be an integer greater than equal to 1.') % display error
elseif (FA * FB) > 0 % Otherwise, if f(a) * f(b) is positive, thereby lacking
                     % evidence that there is a root,
    disp('Error: f(a) * f(b) > 0. Cannot use Bisection Method to determine root. Please input a different [a,b] interval.') % display error.
else % If all of the above conditions are false, we will proceed with the Bisection Method.
    while i <= N % Until we reach the maximum Nth iteration,
        p = a + (b-a)/2; % compute p by finding the midpoint of a and b
        FP = f(p); % compute f(p)
        
        LeftEndPoint(i,1) = a; % Store the ith estimate of endpoint a to the ith entry of the vector
        RightEndPoint(i,1) = b; % Store the ith estimate of endpoint b to the ith entry of the vector
        Approximation(i,1) = p; % Store the ith estimate of the root p to the ith entry of the vector
        FunctionValue(i,1) = FP; % Store the ith estimate of f(p) to the ith entry of the vector
        ErrorBound(i,1) = ((b-a)/2);
        if (FP == 0) || ((b-a)/2 < TOL) % if we've found the root or the size of the ith interval [a,b]
                                        % is less than the input TOL
            
            k = 1; % set k = 1 so the function will be alerted to print the root instead of an error. 
            LeftEndPoint = LeftEndPoint(1:i); % shorten vector so table outputs correctly
            RightEndPoint = RightEndPoint(1:i); % shorten vector
            Approximation = Approximation(1:i); % shorten vector 
            FunctionValue = FunctionValue(1:i); % shorten vector
            ErrorBound = ErrorBound(1:i); % shorten vector
            break   % and then terminate the while loop, and thus the function, since we are done iterating.
        
        else % otherwise, if we have not found the root or reached our tolerance level,
            i = i + 1; % Increment i by one, thus moving onto the next iteration.
            
            if (FA * FP > 0) % If f(a) * f(p) is positive,
                a = p;   % assign the left endpoint to now be p
                FA = FP; % reassign FA to FP.
            else % otherwise, if f(a) * f(p) is negative,
                b = p;  % assign the right endpoint to now be p
            end % close if loop
        
        end % close if loop
    end % close while loop
   
    if k == 1 % if k has been changed from 0 to 1, this means we found a root.
        
        % For table display
        Iteration = [1:i]'; % A vector with entries going from 1 to the last iterate taken. Will be used for tabular output.
        T = table(Iteration, LeftEndPoint, RightEndPoint, Approximation, FunctionValue, ErrorBound); % Make table using appropriate outputs
        disp(T) % Display table
        
        % print the BEST (approximate) root to the 1E-8 as well as 
        % number of iterations taken to find it. This is found by finding
        % the absolute min of the function min table, and then indexing the
        % Approximation (p) column to output this value.
        [~, BestRootIndex] = min(abs(FunctionValue)); % Used to find index of the absolute min in vector FunctionValue, which corresponds to the root
        % print result
        fprintf('Based on the numerical evidence in this table, the best estimate for the root is p = %.8f. TOL reached in %d iterations.\n', Approximation(BestRootIndex,1), i)
        
        % Begin plotting of the endpoints (y-axis) vs. the Iteration number (x-axis)
        scatter(Iteration, LeftEndPoint,'r') % Iteration on x-axis, LeftEndPoint value in red on y-axis
        hold on % hold the figure so both endpoints will be plotted.
        scatter(Iteration, RightEndPoint,'b') % Iteration on x-axis, RightEndpoint in blue on y-axis
        scatter(Iteration, Approximation,'k')
        ylabel('Root Approximation') % label for y axis
        xlabel('Iteration Number') % label for x axis
        legend('Left End Point (a)', 'Right End Point (b)', 'Root Estimate (p)') % Add legend to tell colors apart
        title(func2str(f)) % title for plot is the input function f
        
        % Begin plotting the function y = f(x) on [a-0.25, b+0.25].
        x = linspace(LeftEndPoint(1)-0.25, RightEndPoint(1)+0.25, 101);   % This splits up an interval a little bit larger than 
                                                                % the specified
                                                        % interval into 101
                                                        % points (100 sub-intervals) 
        y = arrayfun(f, x); % compute y = f(x). arrayfun is used to apply a function to each element in an array. 
        figure; % open new figure
        plot(x,y) % plot function
        ylabel('y') % label y-axis
        xlabel('x') % label x-axis
        title(func2str(f)) % title for plot is the input function f
        
    else % if k has not been changed to 1, we did not find the root within
        % the acceptable input tolerance range. Therefore, display the
        % the following error:                   
        fprintf('Bisection Method failed to find an acceptably accurate root within N=%d iterations.\n', N)
    end % close if loop
    
end % close if loop

end % end of function
