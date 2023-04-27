%% Question 25

seed = 1234;
rng(seed)
d = 2; % dimension of the data space
k = 5; % number of klusters
n = 100; % number of data samples
scale = 0.3; % to control separation of the Gaussians klusters
[mu,sigma,w,x] = makedata(d,k,n,scale,false);
[u,X,y] = reparametrize(w,mu,sigma,x);

problem.M = M_MLE3(d,k);
problem.cost = @(point) loglikelyhood(point.u,point.X,y);
problem.egrad = @(point) egrad_l(point.u,point.X,y);
%problem.egrad = @(point) getApproxGradient(problem, point);

checkmanifold(problem.M);
checkgradient(problem);

%% Question 28

clear
disp("Question 28")
% * d=2; k=1;
theta0.w  = 1;
theta0.g{1}.m  = [0;0];
theta0.g{1}.s  = [1,0;0,1];

theta.w  = 1;
theta.g{1}.m  = [1;0];
theta.g{1}.s  = [0.5,0.25;0.25,1];

err1 = Err(theta,theta0)

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

err2 = Err(theta,theta0)

%% Question 29
% code of CGD

%% Question 30

% function [x, cost, info, options] = conjugategradient(problem)

%% Question 31

seed = 1234;
rng(seed)
k = 1;
d = 2;
n = 1000;
scale = 1;
[mu,sigma,w,xx] = makedata(d,k,n,scale,false);
[u,X,y] = reparametrize(w,mu,sigma,xx);

M = M_MLE3(d,k);
problem.M = M;
problem.cost = @(point) loglikelyhood(point.u,point.X,y);
problem.egrad = @(point) egrad_l(point.u,point.X,y);
%problem.egrad = @(point) getApproxGradient(problem, point);

x0 = problem.M.rand();
option.x0 = x0;
option.maxtime = 5;
option.maxiter = Inf;
option.tolgradnorm = 1e-5;
[x, cost, info, option] = RGD(problem, option);
[x1, cost1, info1, options1] = conjugategradient(problem, x0, option);
%
figure;
semilogy([info.iter], [info.gradnorm], '.-');
hold on;
semilogy([info1.iter], [info1.gradnorm], '.-');
legend("Riemanian gradient descent","Conjugated gradient descent")
xlabel('iteration');
ylabel('gradient norm');
title("gradient descents")

figure;
plot([info.iter], [info.cost], '.-');
hold on;
plot([info1.iter], [info1.cost], '.-');
legend("Riemanian gradient descent","Conjugated gradient descent")
xlabel('iteration');
ylabel('cost');
title("gradient descents")
%
figure;
[w,mu,sigma] = deparametrize(x1.u,x1.X);
plot(xx(1,:),xx(2,:),'.','MarkerSize', 8);
hold on;
ell = error_ellipse(sigma{1},mu{1},'conf',0.999)


