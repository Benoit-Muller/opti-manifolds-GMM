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
problem.cost = @(point) loglikelyhood(point.S,point.P,y);
problem.egrad = @(point) egrad_l(point.S,point.P,y);
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

clear
disp("Question 29")
% Generate random problem data.
n = 2;
D = diag(1+rand(n, 1));
[Q, R] = qr(randn(n));
A = Q*D*Q';
% Create the problem structure.
manifold = spherefactory(n);
problem.M = manifold;
% Define the problem cost function and its Euclidean gradient.
problem.cost  = @(x) -x'*(A*x);
problem.egrad = @(x) -2*A*x;
% Solve.
option.maxtime = 10;
option.maxiter = inf;
[x, cost, info, option] = RGD(problem, option);
% Display some statistics.
display(info)
display(cost)
display(-norm(A))
