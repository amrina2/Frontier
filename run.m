%%%%%%
%%%%% File name: run.m 
%%%%% Computing Artifact
%%%%% Author: Amrina Ferdous
%%%%% Purpose: This is the file that researchers need to run only. All
%%%%% other files in this folder are the supporting files due to run.m
%%%%% file.
%%%%%%%%%%%%%%%%%%%%%%%
%
m=100 % Integration grid
W = 100 % Maximum of the grid
wmin=0; % Minimum value of w
wmax=W; % Maximum value of w
w = linspace(0,W,m); % Defining the w 
n=11 % Number of data which could be arbitrary
x = linspace(0,W,n)' ; % Defining the parameter here
%
%Importing the data as .txt file
d = importfile("/Users/amrinaferdous/Desktop/Tools/Git/Frontier/data/data.txt", [1, Inf]);
sigma_d= 0.1 ; % Standard uncertainty of the data
sigma = 5 ; %Prior uncertainty which is changeable and is the reason of having different priori covariance.
theta = 1; % This theta is a parameter of the priori covariance. 
K=1; % Number of iteration that we want. This is also changeable. 
%%%%
% Calling z_hat function i.e. the frontier from the main.m file. For more
% details about z_hat, please see the main.m file.
z_hat=main(w,x,d,sigma_d,sigma,theta,K,@Gfun,@ffun)
%%%%
