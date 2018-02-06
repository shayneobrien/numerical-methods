%% Hermite Interpolation Algorithm

% Name: Shayne O'Brien
% Course: MATH 345 (Dr. Haddad)
% Due Date: Saturday, 10/15/16 by 11:59 pm
% Content: Hermite Interpolation Project Part 1

%% Part 1a
% (a) Write Matlab code to find the Hermite polynomial for the following
% problem. A car traveling along a straight road is clocked at a number of
% points. The data from the observations are given in the following table,
% where time is in seconds, distance is in feet, and speed is in
% feet/second.

% This cell will compute the coefficients of the Hermite Polynomial for the
% data given in data.m, display the coefficients and then display the
% Hermite interpolation polynomial.

% INPUTS: the numbers x0,x1,...,xn; values of f(x0),...,f(xn) and
% f'(x0),...,f'(xn)
% In the case of this project, x = time, f(x) = distance, f'(x) = speed

% OUPUTS: the numbers Q(0,0), Q(1,1), ..., Q(2n+1,2n+1) where
    % H(x) = Q(0,0) + Q(1,1)*(x-x0) + Q(2,2)*(x-x0)^2 +
    % Q(3,3)*(x-x0)^2*(x-x1) + Q(4,4)*(x-x0)^2*(x-x1)^2 + ...
    % + Q(2n+1,2n+1)*(x-x0)^2(x-x1)^2*...*(x-xn-1)^2*(x-xn)
    % In this lab, we will also print the corresponding Hermite polynomial.

format long; % set format to long
data; % load data from data.m file
      % imports three variables: 
      % Time = [0,3,5,8,13]; time in seconds
      % Distance = [0,222,380,625,990]; distance in feet = f(x)
      % Speed = [75,77,80,75,72]; speed in feet per second  = f'(x)

      
% Rename the variables for simplification purposes in coding
x = Time; % let x = Time. We will refer to time values as x from here onwards instead of t.
f = Distance; % let f = f(x) = Distance
df = Speed; % let df = f'(x) = Speed
      
n = length(Time); % set n to the total number of data points we have for x
z = zeros(1,2*n); % let z be a vector of 2*n zeros
Q = zeros(2*n-1,2*n-1); % let Q be a vector of 2*n zeros

for i=1:n % For all x nodes:
    z(2*i-1) = x(i); % Fill vector z with duplicates of each node 
    z(2*i) = x(i); % Example: If data are [0,3,5], z will become [0,0,3,3,5,5].
    Q(2*i-1,1) = f(i); % Computing values for divided differences matrix Q
    Q(2*i,1) = f(i); % Computing values for divided differences matrix Q
    Q(2*i,2) = df(i); % Computing values for divided differences matrix Q
    if i ~= 1 % Special case since for the divided differences matrix the first column is determined by the initial data
        Q(2*i-1,2) = (Q(2*i-1,1)-Q(2*i-2,1))/(z(2*i-1)-z(2*i-2)); % Computing values for divided differences matrix Q
    end % end inner if loop
end % end outer if loop

for i = 2:(2*n-1) % Filling rest of table according to divided differences formula
    for j = 2:i % so that we will go through all possible i,j values
        Q(i+1,j+1) = (Q(i+1,j) - Q(i,j)) ./ (z(i+1)-z(i-j+1)); % Divided differences formula
    end % end inner for loop
end % end outer for loop


% Displaying the coefficients for H(x):
disp('Coefficients of the Hermite interpolating polynomial:')
for i=1:size(Q,1) % For every entry on the diagonal of matrix Q:
    fprintf('Q(%d,%d): %f\n',i,i,Q(i,i)) % print the diagnol entry Q(i,i).
end % end for loop
  

prod = 1; % Initialize product be 1 to start out with
polynomial = Q(1,1); % Q(1,1) = f[x1], which always starts the Hermite interpolation polynomial
syms x; % Let x be a symbol variable so we may display the polynomial
for i = 1:size(Q,1)-1 % until the last entry of the Hermite polynomial:
                      % note: we subtract 1 from end limit of i since we originally set the polynomial as Q(1,1)
    prod = prod*(x-z(i)); % Building polynomial
    newterm = Q(i+1,i+1) * prod; % Multiply polynomial piece by it's coefficient
    polynomial = polynomial + newterm; % Add new piece to existing polynomial
end % end for loop

% print result
fprintf('The Hermite polynomial is therefore computed to be: %s\n\n', char(polynomial))


%% Part 1b
% (b) Use the Hermite polynomial to predict the position of the car at t =
% 12s. (You may calculate this with your code.)

format long; % put format as long
H = matlabFunction(polynomial); % This is a function that converts a sym input to a an function handle
                                % H is now a function that can be passed
                                % inputs. In this case, it is the Hermite
                                % interpolating polynomial H(x).
fprintf('Using the computed hermite polynomial, we estimate the position of the car\nat time t = 12s to be %f feet.\n',H(12))

%% Part 1c
% (c) Use the derivative of the Hermite polynomial to predict whether the
% car ever exceeds 50 miles/hour and predict the maximum speed. Be careful!


% First off, we must change miles/hour into feet/second because our data is
% given in feet/second.

% 50 miles/hour = 50 *(5280 feet/mile)/(60 seconds/min * 60 mins/hour)
%               = 264000 feet / 3600 seconds
%               = 73.333333 feet/second.

% Using our code, we computed the Hermite interpolation polynomial to be:

% H(x) = 75*x 
%         + (2*x^2*(x - 3))/9 
%         - (7*x^2*(x - 3)^2)/225
%         - (29*x^2*(x - 3)^2*(x - 5))/4500 
%         + (163*x^2*(x - 3)^2*(x - 5)^2)/72000 
%         - (2105683025775005*x^2*(x - 3)^2*(x - 5)^2*(x - 8))/2305843009213693952 
%         + (4815588638242799*x^2*(x - 3)^2*(x - 5)^2*(x - 8)^2)/36893488147419103232 
%         - (5968962611917351*x^2*(x - 3)^2*(x - 5)^2*(x - 8)^2*(x - 13))/295147905179352825856

% Computing the derivative of H(x), H'(x):

% H'(x) = 75
%            + (4*x*(x - 3))/9 
%            - (7*x^2*(2*x - 6))/225 
%            - (29*x^2*(x - 3)^2)/4500 
%            - (14*x*(x - 3)^2)/225 
%            + (2*x^2)/9 
%            - (29*x^2*(2*x - 6)*(x - 5))/4500 
%            + (163*x*(x - 3)^2*(x - 5)^2)/36000 
%            + (163*x^2*(2*x - 6)*(x - 5)^2)/72000 
%            + (163*x^2*(2*x - 10)*(x - 3)^2)/72000 
%            - (2105683025775005*x^2*(x - 3)^2*(x - 5)^2)/2305843009213693952 
%            - (29*x*(x - 3)^2*(x - 5))/2250 
%            - (2105683025775005*x*(x - 3)^2*(x - 5)^2*(x - 8))/1152921504606846976 
%            - (2105683025775005*x^2*(2*x - 6)*(x - 5)^2*(x - 8))/2305843009213693952 
%            - (2105683025775005*x^2*(2*x - 10)*(x - 3)^2*(x - 8))/2305843009213693952 
%            + (300974289890175*x*(x - 3)^2*(x - 5)^2*(x - 8)^2)/1152921504606846976 
%            + (300974289890175*x^2*(2*x - 6)*(x - 5)^2*(x - 8)^2)/2305843009213693952 
%            + (300974289890175*x^2*(2*x - 10)*(x - 3)^2*(x - 8)^2)/2305843009213693952 
%            + (300974289890175*x^2*(2*x - 16)*(x - 3)^2*(x - 5)^2)/2305843009213693952 
%            - (5968962611917351*x^2*(x - 3)^2*(x - 5)^2*(x - 8)^2)/295147905179352825856 
%            - (5968962611917351*x*(x - 3)^2*(x - 5)^2*(x - 8)^2*(x - 13))/147573952589676412928 
%            - (5968962611917351*x^2*(2*x - 6)*(x - 5)^2*(x - 8)^2*(x - 13))/295147905179352825856 
%            - (5968962611917351*x^2*(2*x - 10)*(x - 3)^2*(x - 8)^2*(x - 13))/295147905179352825856 
%            - (5968962611917351*x^2*(2*x - 16)*(x - 3)^2*(x - 5)^2*(x - 13))/295147905179352825856

% Next, we sketch the plot of H'(x). The code to do this is located in
% Part 1d. Based  on the plot of H'(x), we see that the maximum speed
% occurs somewhere between time x = 12 and time x = 13. To find the exact
% point, we will first compute the derivative of H'(x), H''(x), and then
% use the Bisection Method to find the critical point of H''(x) which will
% correspond to the maximum of H'(x) = speed. 

% H''(x) is computed by Matlab to be: 

% H''(x) =  (4*x*(x-3))/9-(7*x^2*(2*x-6))/225
%             - (29*x^2*(x-3)^2)/4500
%             - (14*x*(x-3)^2)/225
%             + (2*x^2)/9
%             - (29*x^2*(2*x-6)*(x-5))/4500
%             + (163*x*(x-3)^2*(x-5)^2)/36000
%             + (163*x^2*(2*x-6)*(x-5)^2)/72000
%             + (163*x^2*(2*x-10)*(x-3)^2)/72000
%             - (2105683025775005*x^2*(x-3)^2*(x-5)^2)/2305843009213693952
%             - (29*x*(x-3)^2*(x-5))/2250
%             - (2105683025775005*x*(x-3)^2*(x-5)^2*(x-8))/1152921504606846976
%             - (2105683025775005*x^2*(2*x-6)*(x-5)^2*(x-8))/2305843009213693952
%             - (2105683025775005*x^2*(2*x-10)*(x-3)^2*(x-8))/2305843009213693952
%             + (300974289890175*x*(x-3)^2*(x-5)^2*(x-8)^2)/1152921504606846976
%             + (300974289890175*x^2*(2*x-6)*(x-5)^2*(x-8)^2)/2305843009213693952
%             + (300974289890175*x^2*(2*x-10)*(x-3)^2*(x-8)^2)/2305843009213693952
%             + (300974289890175*x^2*(2*x-16)*(x-3)^2*(x-5)^2)/2305843009213693952
%             - (5968962611917351*x^2*(x-3)^2*(x-5)^2*(x-8)^2)/295147905179352825856
%             - (5968962611917351*x*(x-3)^2*(x-5)^2*(x-8)^2*(x-13))/147573952589676412928
%             - (5968962611917351*x^2*(2*x-6)*(x-5)^2*(x-8)^2*(x-13))/295147905179352825856
%             - (5968962611917351*x^2*(2*x-10)*(x-3)^2*(x-8)^2*(x-13))/295147905179352825856
%             - (5968962611917351*x^2*(2*x-16)*(x-3)^2*(x-5)^2*(x-13))/295147905179352825856
%             + 75

% Letting Hdoubleprime be the function handle for H''(x):
Hdoubleprime = @(x)((4*x)/3 - (28*x*(2*x - 6))/225 - (29*x*(x - 3)^2)/1125 - (29*x^2*(x - 5))/2250 - (14*(x - 3)^2)/225 - (29*x^2*(2*x - 6))/2250 - (29*(x - 3)^2*(x - 5))/2250 + (163*x^2*(x - 3)^2)/36000 + (163*x^2*(x - 5)^2)/36000 + (163*(x - 3)^2*(x - 5)^2)/36000 - (14*x^2)/225 + (300974289890175*(x - 3)^2*(x - 5)^2*(x - 8)^2)/1152921504606846976 + (163*x*(2*x - 6)*(x - 5)^2)/18000 + (163*x*(2*x - 10)*(x - 3)^2)/18000 - (2105683025775005*x*(x - 3)^2*(x - 5)^2)/576460752303423488 - (2105683025775005*x^2*(x - 3)^2*(x - 8))/1152921504606846976 - (2105683025775005*x^2*(x - 5)^2*(x - 8))/1152921504606846976 + (163*x^2*(2*x - 6)*(2*x - 10))/36000 - (2105683025775005*x^2*(2*x - 6)*(x - 5)^2)/1152921504606846976 - (2105683025775005*x^2*(2*x - 10)*(x - 3)^2)/1152921504606846976 - (2105683025775005*(x - 3)^2*(x - 5)^2*(x - 8))/1152921504606846976 + (300974289890175*x^2*(x - 3)^2*(x - 5)^2)/1152921504606846976 + (300974289890175*x^2*(x - 3)^2*(x - 8)^2)/1152921504606846976 + (300974289890175*x^2*(x - 5)^2*(x - 8)^2)/1152921504606846976 - (29*x*(2*x - 6)*(x - 5))/1125 - (2105683025775005*x*(2*x - 6)*(x - 5)^2*(x - 8))/576460752303423488 - (2105683025775005*x*(2*x - 10)*(x - 3)^2*(x - 8))/576460752303423488 - (2105683025775005*x^2*(2*x - 6)*(2*x - 10)*(x - 8))/1152921504606846976 + (300974289890175*x*(2*x - 6)*(x - 5)^2*(x - 8)^2)/576460752303423488 + (300974289890175*x*(2*x - 10)*(x - 3)^2*(x - 8)^2)/576460752303423488 + (300974289890175*x*(2*x - 16)*(x - 3)^2*(x - 5)^2)/576460752303423488 - (5968962611917351*x*(x - 3)^2*(x - 5)^2*(x - 8)^2)/73786976294838206464 - (5968962611917351*x^2*(x - 3)^2*(x - 5)^2*(x - 13))/147573952589676412928 - (5968962611917351*x^2*(x - 3)^2*(x - 8)^2*(x - 13))/147573952589676412928 - (5968962611917351*x^2*(x - 5)^2*(x - 8)^2*(x - 13))/147573952589676412928 + (300974289890175*x^2*(2*x - 6)*(2*x - 10)*(x - 8)^2)/1152921504606846976 + (300974289890175*x^2*(2*x - 6)*(2*x - 16)*(x - 5)^2)/1152921504606846976 + (300974289890175*x^2*(2*x - 10)*(2*x - 16)*(x - 3)^2)/1152921504606846976 - (5968962611917351*x^2*(2*x - 6)*(x - 5)^2*(x - 8)^2)/147573952589676412928 - (5968962611917351*x^2*(2*x - 10)*(x - 3)^2*(x - 8)^2)/147573952589676412928 - (5968962611917351*x^2*(2*x - 16)*(x - 3)^2*(x - 5)^2)/147573952589676412928 - (5968962611917351*(x - 3)^2*(x - 5)^2*(x - 8)^2*(x - 13))/147573952589676412928 - (5968962611917351*x*(2*x - 6)*(x - 5)^2*(x - 8)^2*(x - 13))/73786976294838206464 - (5968962611917351*x*(2*x - 10)*(x - 3)^2*(x - 8)^2*(x - 13))/73786976294838206464 - (5968962611917351*x*(2*x - 16)*(x - 3)^2*(x - 5)^2*(x - 13))/73786976294838206464 - (5968962611917351*x^2*(2*x - 6)*(2*x - 10)*(x - 8)^2*(x - 13))/147573952589676412928 - (5968962611917351*x^2*(2*x - 6)*(2*x - 16)*(x - 5)^2*(x - 13))/147573952589676412928 - (5968962611917351*x^2*(2*x - 10)*(2*x - 16)*(x - 3)^2*(x - 13))/147573952589676412928 - 4/3);

% Running the Bisection Method by using the code we see that there is a
% critical point at the value x = 12.37184029
bisect(Hdoubleprime, 12,13, 1e-13,50)

% Therefore, a maximum occurs at x = 12.37184029. Plugging this into our
% H'(x) and print out the result:
MaxSpeed = Hprime(12.37184029);
fprintf('The maximum speed reached according to our Hermite Polynomial is %f feet/second.\n', MaxSpeed)

% Since the maximum value of the function H'(x) is 119.417339 and H'(x) 
% corresponds to the speed of the car, we thus conclude that the car does 
% in fact exceed 50 mph during the timeframe x = [0,13] since 119.417339 ft/sec > 73.33333 ft/sec = 50 mph.

% We can verify this graphically, which we will do in Part 1d.

%% Part 1d
% (d) Plot the function you get by subdividing the interval where it is
% defined into 100 parts, evaluating the resulting Hermite polynomial at
% these points, and plotting them.

% Plotting the derivative using a function handle and the line y = 50:

clf; % clear any existing plots
% set Hprime as a function handle for H'(x)
Hprime = @(x)((4*x*(x-3))/9-(7*x^2*(2*x-6))/225-(29*x^2*(x-3)^2)/4500-(14*x*(x-3)^2)/225+(2*x^2)/9-(29*x^2*(2*x-6)*(x-5))/4500+(163*x*(x-3)^2*(x-5)^2)/36000+(163*x^2*(2*x-6)*(x-5)^2)/72000+(163*x^2*(2*x-10)*(x-3)^2)/72000-(2105683025775005*x^2*(x-3)^2*(x-5)^2)/2305843009213693952-(29*x*(x-3)^2*(x-5))/2250-(2105683025775005*x*(x-3)^2*(x-5)^2*(x-8))/1152921504606846976-(2105683025775005*x^2*(2*x-6)*(x-5)^2*(x-8))/2305843009213693952-(2105683025775005*x^2*(2*x-10)*(x-3)^2*(x-8))/2305843009213693952+(300974289890175*x*(x-3)^2*(x-5)^2*(x-8)^2)/1152921504606846976+(300974289890175*x^2*(2*x-6)*(x-5)^2*(x-8)^2)/2305843009213693952+(300974289890175*x^2*(2*x-10)*(x-3)^2*(x-8)^2)/2305843009213693952+(300974289890175*x^2*(2*x-16)*(x-3)^2*(x-5)^2)/2305843009213693952-(5968962611917351*x^2*(x-3)^2*(x-5)^2*(x-8)^2)/295147905179352825856-(5968962611917351*x*(x-3)^2*(x-5)^2*(x-8)^2*(x-13))/147573952589676412928-(5968962611917351*x^2*(2*x-6)*(x-5)^2*(x-8)^2*(x-13))/295147905179352825856-(5968962611917351*x^2*(2*x-10)*(x-3)^2*(x-8)^2*(x-13))/295147905179352825856-(5968962611917351*x^2*(2*x-16)*(x-3)^2*(x-5)^2*(x-13))/295147905179352825856+75);

% Begin computing results:
x = linspace(0,13,100); % make a vector of 100 evenly spaced points between 0 and 13
y = zeros(1,100); % let y be a 1x100 vector of zeros
for i = 1:length(x) % For all values of x:
    y(i) = Hprime(x(i)); % evalute H'(x) between [0,13], our defined interval and assign them to the ith entry of vector y
end % end for loop

% Begin plotting results:
plot(x,y,'r') % plot H'(x)
hold on; % hold the figure so we can add additional parts to it
yline = 73.333333*ones(1,100); % let yline be the vector such that y = 50 mph
plot(x,yline,'k') % plot y = 50mph
xlim([0,13]) % adjust x-axis so we only see the defined interval of the problem
xlabel('Time in seconds') % label x-axis
ylabel('Speed in feet/second') % label y-axis
title('Plotting Hermite Polynomial Derivative vs. y = 50 mph = 73.3333 ft/s') % give graph a title

% red line is H'(x)
% black line is line y = 50 mph = 73.33333 ft/s   

% Based on this graph, we can also see that the car does exceed 50 MPH
% during the interval [0,13]. This is a graphical verification of Part 1d.
