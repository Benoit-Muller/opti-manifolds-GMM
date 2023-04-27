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

function X = phi(mu,sigma)
    X = [sigma + mu*mu' , mu;
         mu'            , 1  ];
end