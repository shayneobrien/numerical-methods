%% Principle Component Analysis: Using Covariance Matrix
% Numerical Analysis - Final Project
% David Clarkson, Stephanie Allen, Shayne O'Brien
% Final Project: Due 12/11/16

%% This Program finds axes of maximal variance on a multivariate dataset and projects data points onto them.
% Principle Component analysis is a technique used in exploratory data
% analysis, in which a set of observations of potentially correlated variables
% is transformed into a set of linearly uncorrelated variables, the
% 'principle components.' This is achieved by finding the eigenvectors and
% eigenvalues of the covariance matrix and projecting the data points onto
% the eigenvectors, such that the transformed set of data points is linearly
% uncorrelated (their covariance matrix is diagonal, with covariance terms
% being zero). The principle components are organized by greatest to least
% variance. For highly correlated data points, the first few principle
% components often contain most of the information about each observation,
% and higher order principle components may be left out without much
% information loss, compressing the data. This code presents the user with
% the principle component eigenvalues and allows them to choose how many
% components to use, then plots the first principle component against the
% second principle component.

% NOTE: This method is less efficient than the method using Singular Value
% Decomposition, which does not require the covariance matrix to be
% calculated.

% Main sources: Numerical analysis textbook, Leon textbook (Linear Algebra
% with Applications), Various university power points (which will be cited
% in the report)
%% Data Loading
%Formatting, data loading, space clearing
format longg
clear All
data_NA_Project

%Stores dimensions of data matrix. Note that the rows are observations, and
%the columns are features
[m,n] = size(A);
%Calculates the mean value of each data feature
x_bars = mean(A,1); 

%% Constructing Covariance Matrix
%Preallocates space for X
X = zeros(m,n);

% Creates Matrix X, which contains each (x_i-xbar) and is used to construct
% the covariance matrix
for i = 1:n
    X(1:m,i) = A(1:m,i) - x_bars(i); 
end

% Constructs the covariance matrix. Here m-1 is the number of degrees of
% freedom
COV = (1/(m-1))*(X' * X);

%% Eigenvalues
%Finds the eigenvalues and vectors of the Covariance matrix
[V,D] = eig(COV);

%Sorts the eigenvalues of the covariance matrix by magnitude and stores the
%indices of the associated eigenvectors
[sorted_eig,indices] = sort(diag(D),'descend'); %getting indecies pulled from

%% User Display
sorted_eig %display to the user

%Allows user to reduce the dimensions of the data by only using significant Principle
%Components
number_dimensions = input('Based on the sorted eigenvalues, how many eigenvectors would you like to use for PCA? ');

%preallocates space for feature
component = zeros(n,number_dimensions);

%creates sorted principle component matrix with number of eigenvectors user
%requested
for i = 1:number_dimensions
   component(1:n,i) = V(1:n,indices(i));  
end

%Projects observations onto principle components and transforms
%observations onto Principle Component axes.
Projected_Points=A*component;

%Plots first two components against eachother 
figure
plot(Projected_Points(:,1),Projected_Points(:,2),'Ob'),xlabel('First Principle Component'),
ylabel('Second Principle Component'),title('PCA: First Two Components')
%% Error Analysis
%Generates compressed data matrix using user selected principle components
Compressed=Projected_Points*component';
%Calculates absolute and relative error using the 2 norm
absErr=norm(Compressed-A,2);
relErr=norm(Compressed-A,2)/norm(A,2);
%Displays absolute and relative error to user using two norm
fprintf('\nThe absolute error between the compressed data using %d principle components\nand the original data using the 2 norm is %d\n',number_dimensions,absErr);
fprintf('\nThe relative error between the compressed data using %d principle components\n and the original data using the 2 norm is %d\n',number_dimensions,relErr);