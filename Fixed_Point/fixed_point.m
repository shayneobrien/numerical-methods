%% Fixed-point method

% Name: Shayne O'Brien
% Course: MATH 345 (Dr. Haddad)
% Date: 9/14/16
% Content: Fixed-point method

%% Fixed Point Iteration Algorithm

%A fixed point for a function is a number at which the value of the 
% function does not change when the function is applied. I.e., The number
% p is a fixed point for a given function g if g(p) = p. This algorithm
% will try to solve for p.

%Inputs: function g(x), maximum number of iterations N, tolerance
% threshold TOL.

%Output: Fixed point, or where g(x) = x. A table is output that
% shows iterate number, approximation to the fixed-point for the ith
% iterate, the absolute Cauchy error (distance between consecutive
% iterates), and the relative Cauchy error. A plot will also be generated
% where iterations are mapped against the fixed point estimate at that
% iteration.
 
%%
format long % put format on long for precision and outputs.
close all; clear all; % close all plots, clear all variables.

%% Problem prompt:

% An object of mass 100 kg is release from rest 1000m above the ground and
% allowed to fall under the influence of gravity. Assuming that the force
% due to air resistance is proportional to the velocity of the object with
% proportionality constant k = 10kg/sec, determine when the object will
% strike the ground. The object will strike the grounda t the time T
% determined by the equation:
%              1000 = 98.1*T - 981*(1 - exp(-T/10))
% For the Matlab section of this lab, I will replace T with x.

%% Problem contents for part (b):
%(b) Apply the Fixed-Point Algorithm to approximate the solution of this
%problem with absolute error less than 10^-12 and an initial gues of 10. 
% i. Diary the session and output the results to the Matlab Command Window.
% ii. Provide a table, including the current iterate number the new
% approximation to the fixed-point, Pn, the absolute Cauchy error and the
% relative Cauchy error.
% iii. Plot on one graph (and label) the iteration number vs. the
% iteration, Pn.
% iv. Type in at the command prompt hold on (this will allow you to add to
% this graph later.

%% INPUT FUNCTION g(x) HERE: %%
%Trying to solve for g(x) = x. Manipulate the given function f(x) as needed to
% get x by itself on one side. The other side of the equality is your g(x).
% There are many potential equations for g(x), all we need is one that
% converges (and, preferrably, quickly).

g = @(x)( (1000 + 981*(1-exp(-x/10)))/98.1 ); % Input function

%% Start Fixed-Point Algorithm contents

N = input('Input maximum number of iterations: '); % Ask for max number of iterations
pvec = zeros(N,1); % pre-allocate vector to store each p approximation
pvec(1) = input('Input initial approximation: '); % first entry will be the input initial approximation
TOL = input('Input tolerance: '); % Input tolerance, which we will check against |(pn- p(n-1)) / pn|

AbsCauchyError = zeros(N,1); % pre-allocate vector to store absolute Cauchy error for speed
RelCauchyError = zeros(N,1); % pre-allocate vector to store relative Cauchy error for speed
IterationVec = zeros(N,1); % pre-allocate vector to store iterations for speed

RelCauchyError(1) = NaN; % first entry is NaN because error cannot be computed in the first iteration.
AbsCauchyError(1) = NaN; % first entry is NaN because error cannot be computed in the first iteration.
k = 0; % k is a constant that we will use to simplify the decision making process for our output.
       % if k remains 0 after N iterations, we have not found a fixed point
       % within acceptable TOL. If k = 1, we have and will therefore
       % display all relevant results.

i = 1; % set the beginning iteration

if ((N < 1) || (mod(N,1) ~= 0)) % If the max number of iterations N has been set below 1 or N is not an integer:
    disp('Error: Maximum number of iterations must be an integer greater than equal to 1.') % display error.
elseif TOL <= 0 % Otherwise, if the tolerance is negative:
    disp('Error: TOL must be positive.') % display error.
else % otherwise since both criteria are met, proceed with computation.
    while i <= N % Until we reach the maximum number of iterations:
        p = g(pvec(i)); % compute p-sub-i
        IterationVec(i) = i; % update Iteration vector for i-th entry.
        RelCauchyError(i+1) = abs(p - pvec(i)); % update relative Cauchy error vector. index is i+1 because error cannot be computed in the first iteration.
        AbsCauchyError(i+1) = abs((p - pvec(i))/p); % update absolute Cauchy error vector. index is i+1 because error cannot be computed in the first iteration.
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
    if k == 1 % if we found a fixed point that meets the tolerance criteria within N iterations:
        % displaying results matrix
        ResultsMatrix = [IterationVec, pvec, AbsCauchyError, RelCauchyError]; % Store results in a matrix
        fprintf('    Iteration Number\t      Pn\t    |(Pn - Pn-1)/Pn|\t  |Pn - Pn-1| \n') % Print colum nnames
        disp(ResultsMatrix) % print column contents
        fprintf('The approximate fixed point is p = %.8f, reached in %d iterations.\n', pvec(i),i) % display root and # iterations taken
        
        % plotting results
        plot(IterationVec, pvec) % plot iteration vs. root approximation
        xlabel('Iteration number N') % label x-axis
        ylabel('Fixed-point approximation, Pn') % label y-axis
        title(func2str(g)) % title for plot is the input function g
    else % Otherwise, since a fixed point was not found within a reasonable TOL level:
        fprintf('Fixed Point Method failed to find an acceptably accurate root within N=%d iterations.\n', N) % display error
    end % end if loop
end % end while loop

%% Problem contents for part (c):
%(c) Rewrite this as a root-finding problem for some function f(x) and
% repeat part b using your code for the Bisection Method on the interval,
% [10,30]. Plot and label your approximations for the root on the same graph
% defined in part b. (Your use of hold on in part (b) should allow you to
% output to the same graph if you do not open up a new figure.

%% Comparing Fixed Point Method to the Bisection Method

hold on; % so we can compare the FP Method to the Bisect Method by graphing both on the same graph.
f = @(x)( (98.1*x - 981*(1-exp(-x/10)) - 1000) ); % rewrite function for foot-finding problem using Bisection Method by subtracting 1000 from both sides
                                                  % of the equation to get 0 = 98.1*x - 981(1-exp(-x/10)). We will solve for x
bisect( f, 10, 30, 1e-14, 100); % run Bisection Method, where:
                                       % f = (98.1*x - 981*(1-exp(-x/10)) - 1000)
                                       % [a,b] = [10,30]
                                       % TOL = 1e-14
                                       % max iterations N = 100
