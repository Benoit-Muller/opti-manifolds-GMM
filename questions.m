%% Question 25
disp("––– Question 25 –––")

seed = 1234;
rng(seed)
d = 2; % dimension of the data space
k = 5; % number of klusters
n = 100; % number of data samples
scale = 0.3; % to control separation of the Gaussians klusters

display=false
[mu,sigma,w,xx,n] = makedata(d,k,n,scale,diplay);
[u,X,y,Theta] = reparametrize(w,mu,sigma,xx)

problem.M = M_MLE3(d,k);
problem.cost = @(point) loglikelyhood(point.u,point.X,y);
problem.egrad = @(point) egrad_l(point.u,point.X,y);
%problem.egrad = @(point) getApproxGradient(problem, point);

checkmanifold(problem.M);
checkgradient(problem);

%% Question 28
disp("––– Question 28 –––")
clear
disp("Computation of the total variation distance ")
% * d=2; k=1;
theta0.w  = 1;
theta0.g{1}.m  = [0;0];
theta0.g{1}.s  = [1,0;0,1];

theta.w  = 1;
theta.g{1}.m  = [1;0];
theta.g{1}.s  = [0.5,0.25;0.25,1];

disp("First example")
disp(   "true     : 0.4467")
fprintf("computed : %f\n",Err(theta,theta0))

% * d=2; k=2;
theta0.w  = [0.5;0.5];
theta0.g{1}.m  = [0;0];
theta0.g{2}.m  = [-3;1];
theta0.g{1}.s  = [1,0;0,1];
theta0.g{2}.s  = [0.5,-0.25;-0.25,1];

theta.w  = [0.1;0.9];
theta.g{1}.m  = [1;0];
theta.g{2}.m  = [0;10];
theta.g{1}.s  = [0.5,0.25;0.25,1];
theta.g{2}.s  = [4,0;0,3];

disp("Second example")
disp(   "true     : 1.1228")
fprintf("computed : %f\n",Err(theta,theta0))

%% Question 29
% code of CGD

%% Question 30

% function [x, cost, info, options] = conjugategradient(problem)

%% Question 31
disp("––– Question 31 –––")

clear;
close all;

seed=1234;
rng(seed);
k = 1;
d = 2;
n = 1000;
scale = 1;
question31bis(d,k,n,scale,incldude_f);



