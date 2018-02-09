%% Singular Value Decomposition Final Project

% Name: Stephanie Allen, David Clarkson, Shayne O'Brien
% Course: MATH 345 (Dr. Haddad)
% Due Date: Sunday, 12/11/16 by 12:00 pm
% Content: SVD Decomposition for image compression
% Main sources: Numerical Analysis textbook, MIT Gram Schmidt Slide Show,
% MathWorks Documentation, and Noble and Daniel section

%% SVD Decomposition

% In this script, we use SVD as an application for image compression.
% Input: An image or data for an image
% Output: The same image, but with the size of it reduced.

clear all % clear all variables 
format longg % set format to long without scientific notation

Picture = imread('Tucker_1.jpg'); % insert image here. try to keep it small,
                                  % the code takes a while to run.
pic = double(Picture); % convert data content into double from uint
[m, n, ~] = size(Picture); % store the sizes of the matrix input

% If m < n, SVD cannot be successfully carried out according to our textbook. To handle for this,
% we flip the input image here so m < n is true and then flip it back right before plotting.

if m < n % If m is less than n:
    % flip the picture so m is now greater than n.
    A(:,:,1) = pic(:,:,1)';
    A(:,:,2) = pic(:,:,2)';
    A(:,:,3) = pic(:,:,3)';
    [m, n, ~] = size(A); % store the new sizes of the matrix input
else % Otherwise:
    A = pic; % store it normally
end % end if loop


Ustore = zeros(size(A,1),size(A,1),size(A,3)); % preallocate for speed
V_tstore = zeros(size(A,2),size(A,2),size(A,3)); % preallocate for speed
Sstore = zeros(size(A)); % preallocate for speed
Vstore = V_tstore; % preallocate for speed (V ~= V_t bc of slight roundoff errors)
singstore = zeros(size(A,2),1,size(A,3)); % preallocate for speed

for z = 1:3 % For all three dimensions since images are 3D (red, green, blue values):

    %%%%%%%%%%%%%%%%%
    %%% Finding S %%%
    %%%%%%%%%%%%%%%%%

    [eigenvalues_1,indices] = sort(eig(A(:,:,z)' * A(:,:,z)),'descend'); %finding the eigenvalues and the ordering them
    singular_values = sqrt(eigenvalues_1); % for finding the singular values

    S = zeros(size(A,1),size(A,2)); % preallocate S for speed
    for i = 1:m % For entries 1 to m:
       for j = 1:n % For entries 1 to n:
           if (i == j) % If i is equal to j:
               S(i,j) = singular_values(i); % store i-th singular value in S(i,j). This creates the S matrix
           end % end if loop
       end % end for loop
    end % end for loop
                                       
    %%%%%%%%%%%%%%%%%%%
    %%% Finding V_t %%%
    %%%%%%%%%%%%%%%%%%%
                                       
    % [V,D] = eig(A(:,:,z)) returns diagonal matrix D of eigenvalues and matrix V whose
    % columns are the corresponding right eigenvectors, so that A(:,:,z)*V = V*D
    % Matlab Documentation
                                       
    [V,D1] = eig((A(:,:,z))' * (A(:,:,z))); %V has the eigenvectors 
    V_original = V; % placeholder variable

    for i = 1:length(indices) % For however many eigenvalues there are:
       V(:,indices(i)) = V_original(:,i); % order the eigenvectors according to ordering of eigenvals
    end % end for loop

    V_t = V'; % Transpose V and store it (A = U*S*V')
                                       
    %%%%%%%%%%%%%%%%%%
    %%% Finding U %%%%
    %%%%%%%%%%%%%%%%%%
                                       
    % According to our Numerical Analysis textbook, we need to use the first
    % k nonzero singular values and corresponding eigenvectors to get the first
    % k columns of U

    U = zeros(m,m); % pre allocate for speed

    for i = 1:length(singular_values) % For however many singular values there are:
        U(:,i) = (1/S(i,i)) * A(:,:,z) * V(:,i); % compute U
    end % end for loop

    %Need to pick the next set of vectors such that they are orthog to the
    %columns of U thus far; the textbook suggests do this with Gram-Schmidt

    identity = eye(m); %since U is mxm, we need a total of m linearly independent
    %orthogonal vectors for the columns of U

    U_partial = U(:,1:length(singular_values)); % taking the first k columns of U
    k = length(eigenvalues_1); % let k be the number of eigenvalues

    B = [U_partial identity]; % joining the identity and the partially formed U to do Gram-Schmidt

    Q = zeros(m,(m+k)); % Q will serve as the holder of the columns
    R = zeros((m+k),(m+k)); % R is a placeholder variable

    for j = 1:(k+m) %We are going to apply modified Gram Schmidt to the entirety of B
        v = B(:,j); %v will begin as j-th column of B
        for i = 1:(j-1) % For 1 to j-1:
            R(i,j) = (Q(:,i)'*v) / (Q(:,i)' * Q(:,i)); % modified Gram Schmidt improves accuracy
            v = v - R(i,j)*Q(:,i); % subtract the projection
        end % end for loop
        %v is now perpendicular to all of Q
        Q(:,j) = v; % normalize v to be the unit vector qj
    end % end for loop
                                       
                                    

    
    [ww, vv] = size(Q); % store sizes of Q
    TOL_2 = 1e-11; % set a tolerance
    idx = zeros(1,(vv-ww)); % preallocate for speed
    j = 0; % set j as zero. we will use this to keep track of indices we want to remove

    for i = 1:vv % for all columns of Q
        if norm(Q(:,i)) < TOL_2 % if the col norm is smaller than TOL
            j = j+1; % increase j by 1
            idx(j) = i; % store the index
        end % end 
    end % end

    newQ = zeros(ww,ww); % preallocate for speed

    for i = 1:vv % For all columns of newQ
        if any(i == idx) ~= 1 % If i is not one of the eigenvectors we will throw out:
            newQ(:,i) = Q(:,i); % copy it to Q
        end % end if loop
    end % end for loop

    U = newQ; % we now have U

    [m3,n3] = size(U); % store sizes of U

    %%%%Normalizing%%%%%%%
    for i = 1:(m3) % For all the rows of U:
       norm_1 = norm(U(:,i)); % find the norm
       U(:,i) = (1/norm_1) .* U(:,i); % normalize the row
    end % end for loop
                                       
    Ustore(:,:,z) = U; % store this dimension of U in our storage matrix
    Sstore(:,:,z) = S; % store 
    Vstore(:,:,z) = V; % store
    V_tstore(:,:,z) = V_t; % store
    singstore(:,:,z) = singular_values; % store
    fprintf('Finished dimension %d\n', z) % keep user updated w/ where the code is at, since
                                          % large images can take a while
end % end for loop

%%%%Checking to see how our decomposition compares to A%%%%
U = Ustore; % update U
S = Sstore; % update S
V_t = V_tstore; % update V_t
V = Vstore; % update V (V ~= V_t bc of slight roundoff errors)
singular_values = singstore; % update singular values
test = zeros(size(A)); % set test variable to check error
for i = 1:size(A,3) % For all three dimensions
   test(:,:,i) = U(:,:,i) * S(:,:,i) * V_t(:,:,i); % fill test with SVD vals
end

Rerror = norm(test(:,:,1) - A(:,:,1)); % compute error for Red 
Gerror = norm(test(:,:,2) - A(:,:,2)); % compute error for Green
Berror = norm(test(:,:,3) - A(:,:,3)); % compute error for Blue
fprintf('Red error: %d\nGreen error: %d\nBlue error:%d\n',Rerror,Gerror,Berror') % print errors

%% Image Compression

fprintf('The number of singular values is: %d \n', length(singular_values))
decision = input('Would you like to see them (yes / no)? ','s');
decision = lower(decision); % convert to lowercase
if (strcmp(decision,'yes') == 1) || (strcmp(decision,'y')) % if they said yes
    disp(singular_values) % then display the values
end % end if loop
number_dimensions = input('Based on this, how many would you like to use? ');
A_r = zeros(size(A)); % pre allocate for speed
for j = 1:size(A,3) % For all dimensions of A:
    for i = 1:number_dimensions % For the number that they specified
        A_r(:,:,j) = A_r(:,:,j) + singular_values(i,:,j) * U(:,i,j) * V(:,i,j)'; % compress the image 
    end % end for loop
end % end for loop

figure % open up a figure 
imshow((A_r)/256) % plot the new image! 
                  % We divide by 256 because these are RGB values, but
                  % Matlab needs them to be within [0,1] to plot.