%% Test spd_br, loglikelyhood, egrad_l
clear
disp("Test spd_br, loglikelyhood, egrad_l")
%seed = 1234;
%rng(seed)
d = 2; % dimension of the data space
k = 5; % number of klusters
n = 100; % number of data samples
scale = 0.3; % to control separation of the Gaussians klusters
psi = @sqrt;
[mu,sigma,w,x] = makedata(d,k,n,scale,false);
[u,X,y] = reparametrize(w,mu,sigma,x);

Sphere = spherefactory(k);

spd_br = spd_br_factory(d);
prod_spd_br = powermanifold(spd_br, k);
prod_spd_br = rmfield(prod_spd_br,{'exp'});

problem.M = productmanifold(struct('S',Sphere, 'P',prod_spd_br));

checkmanifold(problem.M);

problem.cost = @(point) loglikelyhood(point.S,point.P,y);
%problem.egrad = @(point) egrad_l(point.S,point.P,y);
problem.egrad = @(point) getApproxGradient(problem, point);

checkgradient(problem);

%% Test RGD

clear
disp("Test RGD")
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
disp(info)
disp(cost)
disp(-norm(A))

%%
% function X = phi(mu,sigma)
%     X = [sigma + mu*mu' , mu;
%          mu'            , 1  ];
% end
% 
% function [u,X,y] = reparametrize(w,mu,sigma,x)
%     [k,~] = size(mu);
%     [~,n] = size(x);
%     u = psi(w);
%     X = cell(1,k);
%     for i=1:k
%         X{i} = phi(mu{i},sigma{i});
%     end
%     y = [x;ones(1,n)];
% end