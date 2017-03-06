%% Hermite Interpolation Part 2 - Snoopy Spline

% Name: Shayne O'Brien
% Course: MATH 345 (Dr. Haddad)
% Due Date: Saturday, 10/22/16 by 11:59 pm
% Content: Part 2: Snoopy spline data file

%%
format short;

% interval [17,27.7]
x = [17, 20, 23, 24, 25, 27, 27.7];
a = [4.5, 7.0, 6.1, 5.6, 5.8, 5.2, 4.1]; % a == f(x)
% b = [3.0, -0.198, -0.609, -0.111, 0.154, -0.401];
% c = [-1.101, 0.035, -0.172, 0.669, -0.403, 0.126];
% d = [-0.126, -0.023, 0.280, -0.357, 0.088, -2.568];
df = [3, -4];

% interval [1,17]
x1 = [1,2,5,6,7,8,10,13,17];
a1 = [3,3.7,3.9,4.2,5.7,6.6,7.1,6.7,4.5];
df1 = [1.0, -0.67];
    
% interval [27.7,30]
x2 = [27.7, 28, 29, 30];
a2 = [4.1, 4.3, 4.1, 3.0];
df2 = [0.33, -1.5];