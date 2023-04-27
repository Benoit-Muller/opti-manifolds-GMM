function [w,mu,sigma] = deparametrize(u,X)
    k = length(u);
    w = u.^2;
    for j=1:k
        mu{j} = X{j}(end,1:end-1);
        sigma{j} = X{j}(1:end-1,1:end-1) - mu{j}*mu{j}';
end