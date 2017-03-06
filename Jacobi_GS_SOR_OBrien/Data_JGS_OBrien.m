%% Algorithms 7.1 (Jacobi), 7.2 (Gauss-Seidel), and 7.3 (SOR)

% Name: Shayne O'Brien
% Course: MATH 345 (Dr. Haddad)
% Due Date: Friday, 11/28/16 by 11:59 pm
% Content: Data for problems % 7.3 #11c, d  (these all involve the same matrix)
                             % 7.3 #16
                             % 7.3 #18 
                             % 7.4 #10 (goes with 7.3 #18)
                             % 7.5 #5c

% The following is data for different problems. Uncomment the data that you want to use.
 
%%%%%%%%%%%%%%%%%%%%% Exercise 7.3.11c %%%%%%%%%%%%%%%%%%%%%

%Use the Gauss-Seidel iterative method to approximate the solution to
%the linear syster with a tolerance of 10^-2 and a maximum of 300
%iterations

% a = [1 0 -1; -.5 1 -.25; 1 -.5 1]; % given matrix
% n = rank(a); % dimension of our matrix
% b = [.2; -1.425; 2]; % given b vector
% XO = [0; 0; 0]; % initial approximation vector

%%%%%%%%%%%%%%%%%%%%% Exercise 7.3.11d %%%%%%%%%%%%%%%%%%%%%

% What happens in part (c) when the system is changed slightly?

% a = [1 0 -2; -.5 1 -.25; 1 -.5 1]; % given matrix
% n = rank(a); % dimension of our matrix
% b = [.2; -1.425; 2]; % given b vector
% XO = [0; 0; 0]; % initial approximation

% Gauss-Seidel no longer converges.

%%%%%%%%%%%%%%%%%%%%% Exercise 7.3.16 %%%%%%%%%%%%%%%%%%%%%

% Suppose that an object can be at any one of n1 equally spaced points x0,
% x1, ..., xn. When an object is at location xi, it is equally likely to
% move to either xi-1 or xi+1 and cannot directly move to any other
% location. Consider the probabilities {Pi} that an object starting at
% location xi will reach the left endpoint before reaching the right endpoint 
% xn. Clearly, P0 = 1 and Pn = 0. Since the object can move to xi only from
% xi-1 or xi+1 and does so with probability 0.5 for each of these
% locations, Pi = 0.5*Pi-1 + 0.5*Pi+1 for each i=1,2,...,n-1. 

% a. Show that the equality in the book holds.

% The equality holds because in the first equation of the system, we have
% x1 - 0.5*x2 = 0. Therefore, x1 = 0.5*x2. 
% Plugging this into the second system of the equation, which is
% -0.5*x1 + x2 - 0.5*x3 = 0, we get that
% -0.5*x1 + x2 - 0.5*x3 
% = -0.5*x2 + x2 - 0.5*x3
% = 0.5*x2 - 0.5*x3 = 0, which is only true when x2 = x3. 
% This train of logic continues all the way through the (n-1)th equation,
% leading to all entries of the b vector in Ax=b being 0 except for the first, 
% so the equality holds.

% b. Solve this system using n = 10,50,100.
% c. Change the probabilities to alpha and 1-alpha for movement to the left
% and right, respetively, and derive the linear system similar to the one
% in part (a).
% d. Repeat part (b) with alpha = 1/3


% For part b, uncomment the n you would like to use

% UNCOMMENT ONE
%n = 10;
%n = 50;
%n = 100;

% UNCOMMENT ONE
%alpha = (1/2); %uncomment this for a probability of (1/2) for movement left
%alpha = (1/3); %uncomment this for a probability of (1/3) for movement left

%for i=1:n-1 % all that these for loops are doing is setting up the matrix for the problem
 %   a(i,i) = 1;
  %  a(i+1,i) = -alpha;
   % a(i,i+1) = -(1-alpha);
%end
%a(n,n) = 1;
%a
%b(1)= alpha;
%for i=2:n
 %   b(i)=0;
%end
%for i=1:n;
 %   XO(i)=0;
%end
%XO = XO';

%%%%%%%%%%%%%%%%%%%%% Exercise 7.3.18, 7.4.10 %%%%%%%%%%%%%%%%%%%%%

%7.3.18a: Explain why the system of equations was reordered.

% The system of equations was reordered so that we have pivots

%7.3.18b: Approximate the solution of the resulting linear system to within
% 0.01 in the l(infinity) norm using as initial approximation the vector
% all of whose entries are 1s with (i) the Jacobi method and (ii) the
% Gauss-Seidel method.


%7.4.10a: Explain why the system of equations was reordered.

% The system of equations was reordered so that we have pivots

%7.4.10b: Approximate the solution of the resulting linear system to within
% 0.01 in the l(infinity) norm using as initial approximation the vector
% all of whose entries are 1s with the SOR method with w = 1.25.

% a = [-1 0 0 sqrt(2)/2 1 0 0 0; 0 -1 0 sqrt(2)/2 0 0 0 0; 0 0 -1 0 0 0 .5 0; 0 0 0 -sqrt(2)/2 0 -1 -.5 0; 0 0 0 0 -1 0 0 1; 0 0 0 0 0 1 0 0; 0 0 0 -sqrt(2)/2 0 0 sqrt(3)/2 0; 0 0 0 0 0 0 -sqrt(3)/2 -1];
% b = [0; 0; 0; 0; 0; 10000; 0; 0];
% XO = [1; 1; 1; 1; 1; 1; 1; 1];
% n = rank(a); % dimension

%%%%%%%%%%%%%%%%%%%%% Exercise 7.5.5c %%%%%%%%%%%%%%%%%%%%%

% (i) Use Gaussian elimination and three-digit rounding arithmetic to
% approximate the solutions to the following linear system. (ii) Then use
% one iteration of iterative refinement to improve the approximation, and
% compare the approximations to the actual solutions.

% (i) We use the following row operations:
% 1.
% A(2,:) = A(2,:) - 14.2/1.19*A(1,:)
% A(4,:) = A(4,:) - 15.3/1.19*A(1,:)
% This gives us
% A = [1.19,2.11,-100,1,1.12;0,-25.3,1200,-12.9,-9.92;0,100,-99.9,1,2.15;0,-27.0,1270,-13.9,-10.2]
% 2.
% A(3,:) = A(3,:) - (100/-25.3)*A(2,:)
% A(4,:) = A(4,:) - (-27/-25.3)*A(2,:)
% This gives us
% A = [1.19,2.11,-100,1,1.12;0,-25.3,1200,-12.9,-9.92;0,0,4640,-50.0,-37.1;0,0,-10.6,-0.133,0.387]
% 3.
% A(4,:) = A(4,:) - (-10.6/4640)*A(3,:)
% This gives us
% A = [1.19,2.11,-100,1,1.12;0,-25.3,1200,-12.9,-9.92;0,0,4640,-50.0,-37.1;0,0,0,-0.247,0.302]

% Backsolving, we get that
% -.247*x4 = 0.302, 
% so x4 = -1.12
% 4640*x3 - 50*x2 = 4640*x3 + 56 = -37.1,
% so x3 = -0.0200
% -25.3*x2 + 1200*x3 - 12.9*x2 = -25.3*x2 + 1200*(-0.0200) - 12.9*(-1.12) = -9.92,
% so x2 = -0.118
% 1.19*x1 + 2.11*x2 - 100*x3 + x4 = 1.19*x1 + 2.11*(-0.118) - 100*(-0.0200) - 1.12 = 1.12
% so x1 = 0.185
% Therefore, our approximation vector is x = [0.185,0.0103,-0.0200,-1.12]'

% for (ii), uncomment this and run the code
% a = [1.19 2.11 -100 1; 14.2 -.122 12.2 -1; 0 100 -99.9 1; 15.3 .11 -13.1 -1];
% b = [1.12; 3.44; 2.15; 4.16];
% XO = [0; 0; 0; 0];
% n = rank(a);
