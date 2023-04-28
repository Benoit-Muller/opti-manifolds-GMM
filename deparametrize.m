function [w,mu,sigma,Theta] = deparametrize(u,X)
    k = length(u);
    w = u.^2;
    mu=cell(k,1);
    sigma = cell(k,1);
    Theta.w = w;
    Theta.g = cell(k,1);
    for j=1:k
        mu{j} = X{j}(end,1:end-1);
        sigma{j} = X{j}(1:end-1,1:end-1) - mu{j}*mu{j}';
        Theta.g{j}.m = mu{j};
        Theta.g{j}.s = sigma{j};
    end
end