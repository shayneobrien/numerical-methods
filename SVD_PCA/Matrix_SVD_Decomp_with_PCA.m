%%%SVD Decomposition with PCA%%%

%%%Authors: Stephanie Allen, David Clarkson, and Shayne O'Brien (with help
%%%from Dr. Haddad)
% Course: MATH 345 (Dr. Haddad)
% Due Date: Sunday, 12/11/16 by 12:00 pm
% Content: SVD Decomposition with PCA
% Main sources: Numerical Analysis textbook, MIT Gram Schmidt Slide Show,
% MathWorks Documentation, Leon textbook (Linear Algebra with
% Applications), and various university power points which will be cited in
% the paper/report.


%%%Set Up: To run the code, make sure that our project folder is your
%%%'directory' in MATLAB and that the data you wish to run the code on is
%%%placed in the data_NA_Project.m script/uncommented in the script.  
%%%Also make sure that the data matrix A is such that m>=n.

%%%Input: Data from the data script (data_NA_Project)

%%%This script take the matrix loaded from the data_NA_Project script
%%%and decomposes the mxn matrix A into the U*S*V^t decomposition, where 
%%%U is mxm, S is mxn, and V^t is nxn.
%%%The script will then perform PCA on the matrix - meaning that it will
%%%use the SVD decomposition to reduce the dimensions.
%%%We assume that the columns of the data matrix A are the variables, and
%%%the rows of the data matrix A are the entries.  Therefore, SVD is
%%%reducing the number of columns 
%%%NOTE: This script will only work for matrices such that m >= n because
%%%we used the material presented in the book, and that is an assumption
%%%they make in Section 9.6.  We didn't have enough time to look at another
%%%algorithm - especially because Section 9.6 does not provide pseudocode.
%%%NOTE: In order for the graph to work, you have to set the number of
%%%dimensions (when asked for the number of dimensions) to 2.

%%%Output: The script outputs the U, S, and V_t matrices of the SVD
%%%Decomposition, the error for the SVD decomposition with regards to the matrix
%%%A, a transformed data matrix, and a graph of the transformed data if the
%%%chosen number of dimensions (when prompted to choose) is set at 2.

clear all
format longg

data_NA_Project  %adding in the data

[m,n] = size(A);


%%%%%Transform the Data by Subtracting the Mean of each Variable from the Variable entries%%%%%%
x_bars = mean(A,1); %need to sum along the columns because those are 
                   %the features (right?)                   

A_original = A;

for i = 1:m
    A(i,1:n) = (A(i,1:n) - x_bars); 
end

%%%%%%%Finding S%%%%%%%%%
%We need order the eigenvalues from absolute largest to smallest, so we need
%to keep tract of this reordering, so that we can reorder the eigenvectors
%accordingly.
[eigenvalues_1,indices] = sort(abs(eig(A' * A)),'descend'); %finding the eigen values 
                                                       %and then ordering them
                                                       %from largest to
                                                       %smallest

singular_values = sqrt(eigenvalues_1); %finding the singular values 

S = zeros(m,n); %S has to be mxn

for i = 1:m %creating the S matrix by putting the singular values on the 'diagonal' of S
   for j = 1:n
       if (i == j) 
           S(i,j) = singular_values(i);      
       else
           S(i,j) = 0;
       end
   end
end

%%%%Finding V_t%%%
%From Matlab Documentation: [V,D] = eig(A) returns diagonal matrix D of 
%eigenvalues and matrix V whose columns are the corresponding right 
%eigenvectors, so that A*V = V*D

[V,D1] = eig(A' * A); %V has the eigenvectors 
V_original = V; %keeping track of the original, unordered V

for i = 1:length(indices)  %we need to order the eigenvectors according 
                           %to the ordering of the eigenvalues
   V(:,indices(i)) = V_original(:,i);
end

V_t = V'; 

%%%%%Finding U%%%%%%
%According to our Numerical Analysis textbook, we need to use the first
%k nonzero singular values and corresponding eigenvectors to get the first
%k columns of U

U = zeros(m,m); %since we know U is mxm 

for i = 1:length(singular_values) 
    U(:,i) = (1/S(i,i)) * A * V(:,i); %formula from page 617 of the textbook
end    

%Need to pick the next set of vectors such that they are orthogonal to the
%columns of U thus far; the textbook suggests do this with Gram-Schmidt

identity = eye(m); %since U is mxm, we need a total of m linearly independent
%orthogonal vectors for the columns of U; these serve as vectors to start
%with for the Gram-Schmidt process

k = length(singular_values);
U_partial = U(:,1:k); %taking the first k columns of U

%The matrix B will allow us to find U by joining the partially formed U
%with 
B = [U_partial identity]; %joining the identity and the partially formed U

Q = zeros(m,(m+k)); %Q will serve as the holder of the orthogonal column vectors
R = zeros((m+k),(m+k)); %R is a placeholder

% Modified Gram-Schmidt Process
for j = 1:(k+m) %We are going to apply Gram Schmidt to the entirety of B
    v = B(:,j); %v will begin as j-th column of B
    for i = 1:(j-1) % for 1 to j-1
        R(i,j) = (Q(:,i)'*v) / (Q(:,i)' * Q(:,i)); %use the updated v in each iteration 
                                                   %according to our
                                                   %sources, "more
                                                   %numerically stable"
        v = v - R(i,j)*Q(:,i); %subtract the projection
                               %Once finish loop
                               %v should now be perpendicular to all of Q
    end 
    Q(:,j) = v; 
end

%%%Removing the columns of 0s from Q (which will become U)%%%%
[m2, n2] = size(Q);

TOL_2 = 1e-10; %NOTE: You might have to adjust this depending upon the 
               %the matrix with which you are working. The code should
               %spit out matrix indices errors if TOL needs to be lowered
             
idx = zeros(1,(n2-m2)); % preallocation for speed
j = 0; % set j as a counter

for i = 1:n2 % for all columns of Q
    if norm(Q(:,i)) < TOL_2 % If the norm of a column is less than the TOL
        j = j+1; % increment j by 1
        idx(j) = i; % if the column is too small, record it
    end
end

newQ = zeros(m2,m2); % for speed

for i = 1:n2 % for all the entries in idx
    if any(i == idx)  ~= 1 % check if i is in idx
        newQ(:,i) = Q(:,i); % copy it
    end
end

U = newQ; % We now have U

[m3,n3] = size(U); % store sizes of U

%%%%Normalizing the U matrix%%%%%%%%
for i = 1:(m3)
   norm_1 = norm(U(:,i)); % find norm
   U(:,i) = (1/norm_1) .* U(:,i);  % normalize
end

%%%%Checking to see how our decomposition compares to A%%%%
test = U * S * V_t; % reconstruct A from U, S, V_t
test_2 = test - A; % compute the error
absolute_error = norm(test_2,2);
relative_error = absolute_error / norm(A,2); 
fprintf('The absolute error for the SVD decomp in representing A according to the 2 norm is %d \n',absolute_error) % display error
fprintf('The relative error for the SVD decomp in representing A according to the 2 norm is %d \n', relative_error)

%%
%%%%%%PCA Transformation%%%%%%%
fprintf('The number of singular values is %d \n', length(singular_values))
decision = input('Would you like to see them (yes / no)? ','s'); % ask user if they wanna see them
decision = lower(decision); % convert to lowercase
if (strcmp(decision,'yes') == 1) || (strcmp(decision,'y')) % if they said yes
    disp(singular_values) % then display the values
end % end if loop
number_dimensions = input('Based on this, how many would you like to use? ');

%%%%%PCA Transformation Multiplication - projecting onto a new set of bases%%%%%%%

V_feature = V(:,(1:number_dimensions)); %this reduces the number of variables

A_transformed = A * V_feature;

%%
%%%%Plotting the Transformed Data%%%%%

if number_dimensions == 2 %only plot if number of dimensions is 2
A_transformed_1 = zeros(m,number_dimensions);

A_transformed_2 = zeros(m,number_dimensions);


for i = 1:m %plotting the transformed data to see if we can spot some grouping 
            %We are plotting according to if the county has a specific type
            %of policy
    if (data(i,3) == 0) %if original data set says has broadbased
        A_transformed_1(i,1:2) = A_transformed(i,:);
    else
        A_transformed_2(i,1:2) = A_transformed(i,:);
    end
end

figure
plot(A_transformed_2(:,1),A_transformed_2(:,2),'oc')
hold on
plot(A_transformed_1(:,1),A_transformed_1(:,2),'ob')
title('Transformed Lower Dimensional Data')
xlabel('Principle Component 1')
ylabel('Principle Component 2')
legend('Has Broadbased Categorical Eligibility','Does NOT have Broadbased Categorical Eligibility')

%Finding the means of the 2 classes
mean1 = mean(A_transformed_1,1);

mean2 = mean(A_transformed_2,1);

disp('The means on the two principle components for the counties with Broad Based eligibility are: ')

mean2

disp('The means for the two principle components for the counties WITHOUT Broad based eligibility are: ')

mean1
end

%%
%Error for PCA
absolute_error_1 = norm((A_transformed * V_feature') - A,2);
relative_error_1 = absolute_error_1 / norm(A,2);

fprintf('\nThe absolute error between the compressed data using %d principle components\nand the original data using the 2 norm is %d\n',number_dimensions,absolute_error_1);
fprintf('\nThe relative error between the compressed data using %d principle components\n and the original data using the 2 norm is %d\n',number_dimensions,relative_error_1);



