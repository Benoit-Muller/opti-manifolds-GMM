%% Test spd_br, loglikelyhood, egrad_l
clear
disp("––– Test spd_br, loglikelyhood, egrad_l –––")

%seed = 1234;
%rng(seed)
d = 2; % dimension of the data space
k = 5; % number of klusters
n = 100; % number of data samples
scale = 0.3; % to control separation of the Gaussians klusters

display=false
[mu,sigma,w,xx,n] = makedata(d,k,n,scale,diplay);
[u,X,y,Theta] = reparametrize(w,mu,sigma,xx)

Sphere = spherefactory(k);

spd_br = spd_br_factory(d);
prod_spd_br = powermanifold(spd_br, k);
prod_spd_br = rmfield(prod_spd_br,{'exp'});

problem.M = productmanifold(struct('u',Sphere, 'X',prod_spd_br));

checkmanifold(problem.M);

problem.cost = @(point) loglikelyhood(point.u,point.X,y);
%problem.egrad = @(point) egrad_l(point.u,point.X,y);
problem.egrad = @(point) getApproxGradient(problem, point);

checkgradient(problem);

%% Test (C)RGD

clear
disp("Test RGD")
% Generate random problem data.
n = 2;
D = diag(randn(n, 1));
[Q, ~] = qr(randn(n));
A = Q*D*Q';
% Create the problem structure.
manifold = spherefactory(n);
problem.M = manifold;
% Define the problem cost function and its Euclidean gradient.
problem.cost  = @(xx) -xx'*(A*xx);
problem.egrad = @(xx) -2*A*xx;
% Solve.
x0 = problem.M.rand();
option.x0 = x0;
option.maxtime = Inf;
option.maxiter = Inf;
option.tolgradnorm = 1e-5;
[x, cost, info, option] = RGD(problem, option);
[x1, cost1, info1, options1] = conjugategradient(problem, x0, option);

% Display some statistics.
disp(info(end))
disp(info1(end))
disp(-norm(A))

figure;
semilogy([info.iter], [info.gradnorm], '.-');
hold on;
semilogy([info1.iter], [info1.gradnorm], '.-');
legend("Riemanian gradient descent", "Conjugated gradient descent")
xlabel('iteration');
ylabel('gradient norm');
title("gradient descents")

figure;
plot([info.iter], [info.cost], '.-');
hold on;
plot([info1.iter], [info1.cost], '.-');
legend("Riemanian gradient descent", "Conjugated gradient descent")
xlabel('iteration');
ylabel('cost');
title("gradient descents")

%% Test question31()
disp("––– Test question31() –––")

clear
close all

k = 1;
d = 2;
n = 1000;
scale = 1;
incldude_f = false;
question31bis(d,k,n,scale,incldude_f);





