clear

seed = 1234;
rng(seed)
d = 2; % dimension of the data space
k = 5; % number of klusters
n = 1000; % number of data samples
scale = 0.3; % to control separation of the Gaussians klusters
psi = @sqrt;
[mu,sigma,w,x] = makedata(d,k,n,scale,false);
[u,X,y] = reparametrize(w,mu,sigma,x);

spd_br = spd_br_factory(d);
Sphere = spherefactory(k);

prod_spd_br = powermanifold(spd_br, k);
prod_spd_br = rmfield(prod_spd_br,{'exp'});
M = productmanifold(struct('S',Sphere, 'P',prod_spd_br));

%%
problem.M = M;
problem.cost = @(point) loglikelyhood(point.S,point.P,y);
problem.egrad = @(point) egrad_l(point.S,point.P,y);

checkgradient(problem);

function X = phi(mu,sigma)
    X = [sigma + mu*mu' , mu;
         mu'            , 1  ];
end

function [u,X,y] = reparametrize(w,mu,sigma,x)
    [k,~] = size(mu);
    [~,n] = size(x);
    u = psi(w);
    X = cell(1,k);
    for i=1:k
        X{i} = phi(mu{i},sigma{i});
    end
    y = [x;ones(1,n)];
end