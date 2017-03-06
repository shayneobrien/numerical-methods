%% Algorithms 6.4, 6.5, 6.6

% Name: Shayne O'Brien
% Course: MATH 345 (Dr. Haddad)
% Due Date: Friday, 11/11/16 by 11:59 pm
% Content: Data for 6.5 #7d #8d 
%                   6.6 #7c #8c
%                   6.6 #7d #8d

%%
% The following is data for different problems. Uncomment the data that you
% wish to use

%% Algorithm 6.4 problems
%Exercise 6.5.7d
    %modify the LU factorization algorithm so that it can be used to solve
    %a linear system and then solve the following
a = [2.1756 4.0231 -2.1732 5.1967; -4.0231 6 0 1.1973; -1 -5.2107 1.1111 0; 6.0235 7 0 -4.1561];
n = size(a,1);
b = [17.102; -6.1593; 3.0004; 0];

%Exercise 6.5.8d
    %modify the LU factorization algorithm so that it can be used to solve
    %a linear system and then solve the following linear system
% a = [2.121 -3.460 0 5.217; 0 5.193 -2.197 4.206; 5.132 1.414 3.141 0;-3.111 -1.732 2.718 5.212];
% n = size(a,1);   %dimensions of our matrix
% b = [1.909; 0; -2.101; 6.824];

%% Algorithm 6.5 problems
%Exercise 6.6.7c
    %modify the LDLt factorization algorithm so that it can be used to
    %solve linear systems. Use it to solve the following
% a = [4 1 -1 0; 1 3 -1 0; -1 -1 5 2; 0 0 2 4];
% n = size(a,1);
% b = [7; 8; -4; 6];

%Exercise 6.6.8c
    %modify the LDLt factorization algorithm so that it can be used to
    %solve linear systems. Use it to solve the following
% a = [4 0 2 1; 0 3 -1 1; 2 -1 6 3; 1 1 3 8];
% n = size(a,1);
% b = [-2; 0; 7; -2];

%% Algorithm 6.6 problems
%Exercise 6.6.7d
    %Use the LLt factorization algorithm to solve the following
% a = [6 2 1 -1; 2 4 1 0; 1 1 4 -1; -1 0 -1 3];
% n = size(a,1);
% b = [0; 7; -1; -2];

%Exercise 6.6.8d
    %Use the LLt factorization algorithm to solve the following
% a = [4 1 1 1; 1 3 0 -1; 1 0 2 1; 1 -1 1 4];
% n = size(a,1);
% b = [2; 2; 1; 1];