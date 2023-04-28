function T = make_Theta(w,mu,sigma)
    k = length(mu);
    T = struct("w",w,"g",cell(k,1));
    for j=1:k
        T.g{j}.m = mu{j};
        T.g{j}.s = sigma{j};
    end
end