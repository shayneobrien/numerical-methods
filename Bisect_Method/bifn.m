%% Bisection Method for Exercises 2.1: 5,9,11,13

% Name: Shayne O'Brien
% Course: MATH 345 (Dr. Haddad)
% Date: 9/6/16
% Content: Running Bisection Method for each exercise.

%% Exercise 2.1.9: a
% Sketch the graphs of y = e^x - 2 and y = cos(e^x - 2)

    close all; clear all; %This will be used throughout this .m file to close old plots
                          % and clear old variables.
    
    % plotting y = e^x - 2
    x = linspace(-3, 3, 100); % generate 100 linearly equally spaced points between -3
                              % and 3 in a vector and assign the vector to variable x.
    y = exp(x) - 2; % create a vector y that is e^x - 2 for each entry in vector x
    plot(x,y) % plot x versus y
    title('Plot of y = e^x - 2') % title for plot
    ylabel('y') % label for y axis
    xlabel('x') % label for x axis

    % plotting y = cos(e^x - 2)
    figure; % open new plot so the one for y = e^x - 2
    hold on % necessary so the plot will actually plot
    x1 = linspace(-3,3,100);% generate 100 linearly equally spaced points between -3
                              % and 3 in a vector and assign the vector to variable x1.
    y1 = cos(exp(x1) - 2); % create a vector y1 that is cos(e^x - 2) for each entry in vector x1
    plot(x1,y1) % plot x versus y
    title('Plot of y = cos(e^x - 2)') % title for plot
    ylabel('y') % label for y axis
    xlabel('x') % label for x axis

    
    
%% Exercise 2.1.9: b
% Use the Bisection Method to find an approximation to within 10^-5 to a
% value in [0.5, 1.5] with e^x - 2 = cos(e^x - 2).

    % Using algebra, we can subtract cos(e^x - 2) from both sides of the
    % equality to get e^x - 2 - cos(e^x - 2) = 0. Having this equality, we
    % can now proceed with the computation:
    close all; clear all; % close plots and clear variables from previous problem.
    
    bisect( @(x)exp(x)-2-cos(exp(x)-2), 0.5, 1.5, 10^-5, 25)
    % Answer: The root is approximately p = 1.00762177.


%% Exercise 2.1.11
% Let f(x) = (x+2)(x+1)x(x-1)^3(x-2). To which zero does the Bisection
% Method converge when applied on the following intervals?

%% Exercise 2.1.11: a
% [-3, 2.5]
    close all; clear all;
    bisect( @(x)(x+2)*(x+1)*x*((x-1)^3)*(x-2), -3, 2.5, 10^-5, 30)
    % Answer: The root is approximately p = 2.00000048.

%% Exercise 2.1.11: b
% [-2.5, 3]
    close all; clear all;
    bisect( @(x)(x+2)*(x+1)*x*((x-1)^3)*(x-2), -2.5, 3, 10^-5, 30)
    % Answer: The root is approximately p = -2.00000048.

%% Exercise 2.1.11: c
% [-1.75, 1.5]
    close all; clear all;
    bisect( @(x)(x+2)*(x+1)*x*((x-1)^3)*(x-2), -1.75, 1.5, 10^-5, 30)
    % Answer: The root is approximately p = -0.99999714.

%% Exercise 2.1.11: d
% [-1.5, 1.75]
    close all; clear all;
    bisect( @(x)(x+2)*(x+1)*x*((x-1)^3)*(x-2), -1.5, 1.75, 10^-5, 30)
    % Answer: The root is approximately p = 0.99999714.

%% Exercise 2.1.13
% Find an approximation to 25^(1/3) to within 10^(-4) using the Bisection
% method.

    % To find an approximation to 25^(1/3), we will consider f(x) = x^3-25 = 0.
    % By adding 25 to both sides and then taking the cube root, we get 
    % x = 25^(1/3). We will consider the interval [2,3]. Then,
    close all; clear all;
    bisect( @(x)x^3 - 25, 2, 3, 10^-5, 30)
    % Answer: The root is approximately p = 2.92401886. Therefore, an
    % approximation to 25^(1/3) is 2.92401886.
