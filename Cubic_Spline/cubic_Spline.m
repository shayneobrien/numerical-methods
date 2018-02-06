%% Hermite Interpolation Part 2 - Snoopy Spline

% Name: Shayne O'Brien
% Course: MATH 345 (Dr. Haddad)
% Due Date: Saturday, 10/22/16 by 11:59 pm
% Content: Algorithm 3.5 code, computation and graphing for the spline

%% Part 2a
% Approximate the upper portion of Snoopy using clamped spine interpolants.
% The curve is drawn on a grid from which the table is constructed. use
% Algorithm 3.5 to construct the three clamped cubic splines. Be sure to
% output the given data, a, b, c, d in a table (b, c, d should have at
% least 4 decimal places of accuracy) for each of the three segments. Use
% format short. You may incorporate the data into an m-file or have it
% input into the program.

format short; % set format to 4 decimal places
data; % load in the data from data.m Matlab file

% This code will construct the cubic spline interpolant S for the function
% f defined at the numbers x1 < x2 < ... < xn, satisfying 
% S'(x1) = f'(x1) and S'(xn) = f'(xn)

% INPUTS: n
%         x1, x2, ..., xn
%         a1 = f(x1), a2 = f(x2), ..., an = f(xn)
%         FPO = f'(x1)
%         FPN = f'(xn)

% OUTPUTS: aj, bj, cj, dj for j = 1, ..., n-1

% Note: S(x) = aj + bj*(x-xj) + cj*(x-xj)^2 + dj*(x-xj)^3 for xj <= x <= x(j+1)

% We will compute the values of the spline for [17,27.7] in this cell.

n = length(x); % let n be the length of the vector x
h = zeros(1,n-1); % let h be a vector of length n-1 filled with zeros
alpha = zeros(1,n); % let alpha be a vector of length n filled with zeros
l = zeros(1,n); % let l be a vector of length n filled with zeros
u = zeros(1,n-1); % let u be a vector of length n-1 filled with zeros
z = zeros(1,n); % let z be a vector of length n filled with zeros
iteration = 1:n; % let iteration be a vector integers 1,2,...,n

b = zeros(1,n-1);
c = zeros(1, n-1);
d = zeros(1, n-1);

for i = 1:n-1 % for all entries except the last:
    h(i) = x(i+1)-x(i); % compute the i-th entry of h
end % end for loop

alpha(1) = 3*(a(2)-a(1))/h(1) - 3*df(1); % compute the first entry of alpha
alpha(n) = 3*df(2)-3*(a(n)-a(n-1))/h(n-1); % compute the last entry of alpha

for i = 2:n-1 % for the second entry through the penultimate entry:
    alpha(i) = (3/h(i))*(a(i+1) - a(i)) - (3/h(i-1))*(a(i)-a(i-1)); % compute the i-th entry of alpha
end % end for loop

l(1) = 2*h(1); % compute first entry of vector l
u(1) = 0.5; % compute first entry of vector u
z(1) = alpha(1) / l(1); % compute first entry of vector z

% the next two for loops are solving the tridiagnol linear system
for i = 2:n-1 % for the second entry through the penultimate entry:
    l(i) = 2*(x(i+1) - x(i-1)) - h(i-1)*u(i-1); % compute the i-th entry of l
    u(i) = h(i)/l(i); % compute the i-th entry of u
    z(i) = (alpha(i) - h(i-1)*z(i-1))/l(i); % compute the i-th entry of z
end % end for loop

l(n) =  h(n-1)*(2 - u(n-1)); % compute the last entry of l
z(n) = (alpha(n) - h(n-1)*z(n-1))/l(n); % compute the last entry of z
c(n) = z(n); % let the last entry of c be equal to the last entry of z

for j = n-1:-1:1 % going backwards from the penultimate entry to the first entry:
    c(j) = z(j) - u(j)*c(j+1); % compute j-th entry of vector c
    b(j) = (a(j+1) - a(j))/h(j) - h(j)*(c(j+1) + 2*c(j))/3; % compute j-th entry of vector b
    d(j) = (c(j+1) - c(j))/(3*h(j)); % compute j-th entry of vector d
end % end for loop

% we will now move on to actually computing the snoopy spline.
linearvec = linspace(17, 27.7, 1000); % make a vector of 1000 equally spaced points between 17 and 27.7
Sx = zeros(1,length(linearvec)); % make a 1x1000 vector to store values of the curve

% Recall that 
% S(x) = aj + bj*(x-xj) + cj*(x-xj)^2 + dj*(x-xj)^3 for xj <= x <= x(j+1)

for i = 1:length(linearvec)
    if (17 <= linearvec(i)) && (linearvec(i) <= 20)
        j = 1;
        Sx(i) = a(j) + b(j)*(linearvec(i)-x(j)) + c(j)*(linearvec(i)-x(j))^2 + d(j)*(linearvec(i)-x(j))^3;
    elseif (20 <= linearvec(i)) && (linearvec(i) <= 23)
        j = 2;
        Sx(i) = a(j) + b(j)*(linearvec(i)-x(j)) + c(j)*(linearvec(i)-x(j))^2 + d(j)*(linearvec(i)-x(j))^3;
    elseif (23 <= linearvec(i)) && (linearvec(i) <= 24)
        j = 3;
        Sx(i) = a(j) + b(j)*(linearvec(i)-x(j)) + c(j)*(linearvec(i)-x(j))^2 + d(j)*(linearvec(i)-x(j))^3;
    elseif (24 <= linearvec(i)) && (linearvec(i) <= 25)
        j = 4;
        Sx(i) = a(j) + b(j)*(linearvec(i)-x(j)) + c(j)*(linearvec(i)-x(j))^2 + d(j)*(linearvec(i)-x(j))^3;
    elseif (25 <= linearvec(i)) && (linearvec(i) <= 27)
        j = 5;
        Sx(i) = a(j) + b(j)*(linearvec(i)-x(j)) + c(j)*(linearvec(i)-x(j))^2 + d(j)*(linearvec(i)-x(j))^3;
    elseif (27 <= linearvec(i)) && (linearvec(i) <= 27.7)
        j = 6;
        Sx(i) = a(j) + b(j)*(linearvec(i)-x(j)) + c(j)*(linearvec(i)-x(j))^2 + d(j)*(linearvec(i)-x(j))^3;
    end
end

% Display table of results
fprintf('Data input into program for interval [%.0f,%.0f]: \n\n', x(1),x(n)) % tell user what data is being used and skip two lines
DataMatrix = [iteration; x; a]'; % store data into DataMatrix variable
fprintf('       i\tx\t f(x)\n') % print column names and skip a line
disp(DataMatrix) % display the data
fprintf('With FPO = df(x1) = %.4f and FPN = df(xn) = %.4f\n', df(1),df(2)) % tell user what derivatives were

fprintf('Results: \n\n') % tell user we are printing the results, skip two lines
ResultsMatrix = [a(1:n-1); b(1:n-1); c(1:n-1); d(1:n-1)]'; % Store results in a matrix for j = 1, ..., n-1
fprintf('     a=f(x)     b\t  c\t    d\t\n') % Print column names and skip a line
disp(ResultsMatrix) % print column contents

%% Note to reader
% NOTE: The code for the next two parts will look extremely similar
% to the code in this cell. The difference between this cell and the next two cells is
% that we are computing the spline for the other two sections of the curve,
% [1,17] and [27.7,30]. All results will be plotted in the final cell of
% this file. I store the values of the spline for the next two intervals in
% Sx1 and Sx2.

%% Part 2a: interval [1,17]

% We will compute the values of the spline for [1,17] in this cell.

x = x1; % let x be x1, the data for the interval [1,17]
a = a1; % let a be a1, the data for the interval [1,17]
df = df1; % let df be df1, the data for the interval [1,17]

n = length(x); % let n be the length of the vector x
h = zeros(1,n-1); % let h be a vector of length n-1 filled with zeros
alpha = zeros(1,n); % let alpha be a vector of length n filled with zeros
l = zeros(1,n); % let l be a vector of length n filled with zeros
u = zeros(1,n-1); % let u be a vector of length n-1 filled with zeros
z = zeros(1,n); % let z be a vector of length n filled with zeros
iteration = 1:n; % let iteration be a vector integers 1,2,...,n

b = zeros(1,n-1);
c = zeros(1, n-1);
d = zeros(1, n-1);

for i = 1:n-1 % for all entries except the last:
    h(i) = x(i+1)-x(i); % compute the i-th entry of h
end % end for loop

alpha(1) = 3*(a(2)-a(1))/h(1) - 3*df(1); % compute the first entry of alpha
alpha(n) = 3*df(2)-3*(a(n)-a(n-1))/h(n-1); % compute the last entry of alpha

for i = 2:n-1 % for the second entry through the penultimate entry:
    alpha(i) = (3/h(i))*(a(i+1) - a(i)) - (3/h(i-1))*(a(i)-a(i-1)); % compute the i-th entry of alpha
end % end for loop

l(1) = 2*h(1); % compute first entry of vector l
u(1) = 0.5; % compute first entry of vector u
z(1) = alpha(1) / l(1); % compute first entry of vector z

% the next two for loops are solving the tridiagnol linear system
for i = 2:n-1 % for the second entry through the penultimate entry:
    l(i) = 2*(x(i+1) - x(i-1)) - h(i-1)*u(i-1); % compute the i-th entry of l
    u(i) = h(i)/l(i); % compute the i-th entry of u
    z(i) = (alpha(i) - h(i-1)*z(i-1))/l(i); % compute the i-th entry of z
end % end for loop

l(n) =  h(n-1)*(2 - u(n-1)); % compute the last entry of l
z(n) = (alpha(n) - h(n-1)*z(n-1))/l(n); % compute the last entry of z
c(n) = z(n); % let the last entry of c be equal to the last entry of z

for j = n-1:-1:1 % going backwards from the penultimate entry to the first entry:
    c(j) = z(j) - u(j)*c(j+1); % compute j-th entry of vector c
    b(j) = (a(j+1) - a(j))/h(j) - h(j)*(c(j+1) + 2*c(j))/3; % compute j-th entry of vector b
    d(j) = (c(j+1) - c(j))/(3*h(j)); % compute j-th entry of vector d
end % end for loop

% we will now move on to actually computing the snoopy spline.
linearvec1 = linspace(1, 17, 1000); % make a vector of 1000 equally spaced points between 17 and 27.7
Sx1 = zeros(1,length(linearvec1)); % make a 1x1000 vector to store values of the curve

% Recall that 
% S(x) = aj + bj*(x-xj) + cj*(x-xj)^2 + dj*(x-xj)^3 for xj <= x <= x(j+1)

for i = 1:length(linearvec1)
    if (1 <= linearvec1(i)) && (linearvec1(i) <= 2)
        j = 1;
        Sx1(i) = a(j) + b(j)*(linearvec1(i)-x(j)) + c(j)*(linearvec1(i)-x(j))^2 + d(j)*(linearvec1(i)-x(j))^3;
    elseif (2 <= linearvec1(i)) && (linearvec1(i) <= 5)
        j = 2;
        Sx1(i) = a(j) + b(j)*(linearvec1(i)-x(j)) + c(j)*(linearvec1(i)-x(j))^2 + d(j)*(linearvec1(i)-x(j))^3;
    elseif (5 <= linearvec1(i)) && (linearvec1(i) <= 6)
        j = 3;
        Sx1(i) = a(j) + b(j)*(linearvec1(i)-x(j)) + c(j)*(linearvec1(i)-x(j))^2 + d(j)*(linearvec1(i)-x(j))^3;
    elseif (6 <= linearvec1(i)) && (linearvec1(i) <= 7)
        j = 4;
        Sx1(i) = a(j) + b(j)*(linearvec1(i)-x(j)) + c(j)*(linearvec1(i)-x(j))^2 + d(j)*(linearvec1(i)-x(j))^3;
    elseif (7 <= linearvec1(i)) && (linearvec1(i) <= 8)
        j = 5;
        Sx1(i) = a(j) + b(j)*(linearvec1(i)-x(j)) + c(j)*(linearvec1(i)-x(j))^2 + d(j)*(linearvec1(i)-x(j))^3;
    elseif (8 <= linearvec1(i)) && (linearvec1(i) <= 10)
        j = 6;
        Sx1(i) = a(j) + b(j)*(linearvec1(i)-x(j)) + c(j)*(linearvec1(i)-x(j))^2 + d(j)*(linearvec1(i)-x(j))^3;
    elseif (10 <= linearvec1(i)) && (linearvec1(i) <= 13)
        j = 7;
        Sx1(i) = a(j) + b(j)*(linearvec1(i)-x(j)) + c(j)*(linearvec1(i)-x(j))^2 + d(j)*(linearvec1(i)-x(j))^3;
    elseif (13 <= linearvec1(i)) && (linearvec1(i) <= 17)
        j = 8;
        Sx1(i) = a(j) + b(j)*(linearvec1(i)-x(j)) + c(j)*(linearvec1(i)-x(j))^2 + d(j)*(linearvec1(i)-x(j))^3;
    end
end

% Display table of results
fprintf('Data input into program for interval [%.0f,%.0f]: \n\n', x(1),x(n)) % tell user what data is being used and skip two lines
DataMatrix = [iteration; x; a]'; % store data into DataMatrix variable
fprintf('       i\tx\t f(x)\n') % print column names and skip a line
disp(DataMatrix) % display the data
fprintf('With FPO = df(x1) = %.4f and FPN = df(xn) = %.4f\n', df(1),df(2)) % tell user what derivatives were

fprintf('Results: \n\n') % tell user we are printing the results, skip two lines
ResultsMatrix = [a(1:n-1); b(1:n-1); c(1:n-1); d(1:n-1)]'; % Store results in a matrix for j = 1, ..., n-1
fprintf('     a=f(x)     b\t  c\t    d\t\n') % Print column names and skip a line
disp(ResultsMatrix) % print column contents

%% Part 2a: interval [27.7,30]

% We will compute the values of the spline for [17,27.7] in this cell.

x = x2; % let x be x2, the data for the interval [17,27.7]
a = a2; % let a be a2, the data for the interval [17,27.7]
df = df2; % let df be df2, the data for the interval [17,27.7]

n = length(x); % let n be the length of the vector x
h = zeros(1,n-1); % let h be a vector of length n-1 filled with zeros
alpha = zeros(1,n); % let alpha be a vector of length n filled with zeros
l = zeros(1,n); % let l be a vector of length n filled with zeros
u = zeros(1,n-1); % let u be a vector of length n-1 filled with zeros
z = zeros(1,n); % let z be a vector of length n filled with zeros
iteration = 1:n; % let iteration be a vector integers 1,2,...,n

b = zeros(1, n-1);
c = zeros(1, n-1);
d = zeros(1, n-1);

for i = 1:n-1 % for all entries except the last:
    h(i) = x(i+1)-x(i); % compute the i-th entry of h
end % end for loop

alpha(1) = 3*(a(2)-a(1))/h(1) - 3*df(1); % compute the first entry of alpha
alpha(n) = 3*df(2)-3*(a(n)-a(n-1))/h(n-1); % compute the last entry of alpha

for i = 2:n-1 % for the second entry through the penultimate entry:
    alpha(i) = (3/h(i))*(a(i+1) - a(i)) - (3/h(i-1))*(a(i)-a(i-1)); % compute the i-th entry of alpha
end % end for loop

l(1) = 2*h(1); % compute first entry of vector l
u(1) = 0.5; % compute first entry of vector u
z(1) = alpha(1) / l(1); % compute first entry of vector z

% the next two for loops are solving the tridiagnol linear system
for i = 2:n-1 % for the second entry through the penultimate entry:
    l(i) = 2*(x(i+1) - x(i-1)) - h(i-1)*u(i-1); % compute the i-th entry of l
    u(i) = h(i)/l(i); % compute the i-th entry of u
    z(i) = (alpha(i) - h(i-1)*z(i-1))/l(i); % compute the i-th entry of z
end % end for loop

l(n) =  h(n-1)*(2 - u(n-1)); % compute the last entry of l
z(n) = (alpha(n) - h(n-1)*z(n-1))/l(n); % compute the last entry of z
c(n) = z(n); % let the last entry of c be equal to the last entry of z

for j = n-1:-1:1 % going backwards from the penultimate entry to the first entry:
    c(j) = z(j) - u(j)*c(j+1); % compute j-th entry of vector c
    b(j) = (a(j+1) - a(j))/h(j) - h(j)*(c(j+1) + 2*c(j))/3; % compute j-th entry of vector b
    d(j) = (c(j+1) - c(j))/(3*h(j)); % compute j-th entry of vector d
end % end for loop

% we will now move on to actually computing the snoopy spline.
linearvec2 = linspace(27.7, 30, 1000); % make a vector of 1000 equally spaced points between 17 and 27.7
Sx2 = zeros(1,length(linearvec2)); % make a 1x1000 vector to store values of the curve

% Recall that 
% S(x) = aj + bj*(x-xj) + cj*(x-xj)^2 + dj*(x-xj)^3 for xj <= x <= x(j+1)

for i = 1:length(linearvec2)
    if (27.7 <= linearvec2(i)) && (linearvec2(i) <= 28)
        j = 1;
        Sx2(i) = a(j) + b(j)*(linearvec2(i)-x(j)) + c(j)*(linearvec2(i)-x(j))^2 + d(j)*(linearvec2(i)-x(j))^3;
    elseif (28 <= linearvec2(i)) && (linearvec2(i) <= 29)
        j = 2;
        Sx2(i) = a(j) + b(j)*(linearvec2(i)-x(j)) + c(j)*(linearvec2(i)-x(j))^2 + d(j)*(linearvec2(i)-x(j))^3;
    elseif (29 <= linearvec2(i)) && (linearvec2(i) <= 30)
        j = 3;
        Sx2(i) = a(j) + b(j)*(linearvec2(i)-x(j)) + c(j)*(linearvec2(i)-x(j))^2 + d(j)*(linearvec2(i)-x(j))^3;
    end
end

% Display table of results
fprintf('Data input into program for interval [%.1f,%.0f]: \n\n', x(1),x(n)) % tell user what data is being used and skip two lines
DataMatrix = [iteration; x; a]'; % store data into DataMatrix variable
fprintf('       i\tx\t f(x)\n') % print column names and skip a line
disp(DataMatrix) % display the data
fprintf('With FPO = df(x1) = %.4f and FPN = df(xn) = %.4f\n', df(1),df(2)) % tell user what derivatives were

fprintf('Results: \n\n') % tell user we are printing the results, skip two lines
ResultsMatrix = [a(1:n-1); b(1:n-1); c(1:n-1); d(1:n-1)]'; % Store results in a matrix for j = 1, ..., n-1
fprintf('     a=f(x)     b\t  c\t    d\t\n') % Print column names and skip a line
disp(ResultsMatrix) % print column contents

Sn%% Part 2b
% Plot the results on a single graph. In order to plot 3 separate graphs on
% one set of axes, use hold on; after each set is plotted. At the end of
% your code, put in a hold off;

clf; % clear any existing plots

plot(linearvec, Sx, 'g') % plot spline for [17,27.7] in green
hold on; % hold on so we can add another plot
plot(linearvec1, Sx1, 'b') % plot spline for [1,17] in blue
hold on;
plot(linearvec2, Sx2, 'r') % plot spline for [27.7,30] in red
title('Snoopy Spline') % title the graph
xlabel('x') % label x-axis
ylabel('y') % label y-axis
xlim([0, 30]) % set x-axis limits
ylim([0,8]) % set y-axis limits
hold off; % turn off hold
