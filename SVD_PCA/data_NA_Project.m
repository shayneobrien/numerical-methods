%%%%data Script: data_NA_Project%%%%%

%%%%Test Data%%%%%%

%A = randi([1,100],120,40);

%%%%%%%PCA Data - All of it!%%%%%%%
% data = csvread('SNAP_Data_2.csv');
% 
% A = [data(:,1:2) data(:,5) data(:,9:10)];


%%%%PCA Data - Some of it%%%%%%%%
%If you want the run time to be shorter, here is a smaller amount of data
data = csvread('SNAP_Data_2.csv');

data = data(221:720,:); %500 data values

A = [data(:,1:2) data(:,5) data(:,10)]; 










