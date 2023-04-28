clc
seed = 9876;
rng(seed)

%% Question 25
fprintf("\n––– Question 25 –––\n")

clear;
close all;

d = 2; % dimension of the data space
k = 5; % number of klusters
n = 100; % number of data samples
scale = 0.3; % to control separation of the Gaussians klusters

problem = problem_MLE3(d,k,n,scale);

checkmanifold(problem.M);
checkgradient(problem);
saveas(gcf,'graphics/q25_checkgradient.pdf')

%% Question 28
fprintf("\n––– Question 28 –––\n")

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
fprintf("\n––– Question 29 –––\n")

problem = problem_MLE3(2,5,100,0.3);
option = struct("maxiter",Inf, "maxtime",60, "tolgradnorm",1e-3, ...
                "verbosity",0,"x0",problem.M.rand());
[x, cost, info, option] = RGD(problem, option);
disp(info(end))

%% Question 30
fprintf("\n––– Question 30 –––\n")

[x, cost, info, option] = conjugategradient(problem);
disp(info(end))

%% Question 31
seed = 9876;
rng(seed)
fprintf("\n––– Question 31 –––\n")
clear;
close all;

k = 1;
d = 2;
n = 1000; %1000
scale = 1;
incldude_f = true;

question31(d,k,n,scale,incldude_f);

%% Question 32
fprintf("\n––– Question 32 –––\n")
clear;
close all;

kk = 2:5;
d = 2;
n = 1000;
scale = 0.3;
incldude_f = true;

for k = kk
    [r,c] = question31(d,k,n,scale,incldude_f);
end

%% Question 33
fprintf("\n––– Question 33 –––\n")
clear;
close all;

k = [2,4];
dd = [5,10,15];
n = 1000;
scale = 0.3;
incldude_f = false;

for d = dd
    [r,c] = question31(d,k,n,scale,incldude_f);
end