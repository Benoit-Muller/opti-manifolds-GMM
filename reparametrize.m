function [u,X,y,Theta] = reparametrize(w,mu,sigma,x)
    [k,~] = size(mu);
    [~,n] = size(x);
    u = psi(w);
    X = cell(1,k);
    for i=1:k
        X{i} = phi(mu{i},sigma{i});
    end
    y = [x;ones(1,n)];
    Theta.u=u;
    Theta.X =X;
    Theta.w=w;
    Theta.m=mu;
    Theta.s=sigma;
end

function X = phi(mu,sigma)
    X = [sigma + mu*mu' , mu;
         mu'            , 1  ];
end