function pp = p(mu,sigma,x)
% gaussian density p_d
    d=length(mu);
    pp = (det(sigma)*(2*pi)^d)^(-0.5) * exp(-0.5*(x-mu)'*(sigma\(x-mu)));
end