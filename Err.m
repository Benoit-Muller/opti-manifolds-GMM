function err = Err(theta, theta0)
    k = length(theta.w);
    permutations = perms(1:k);
    err = inf;
    for i =1:factorial(k)
        err = min([err, Errp(theta, theta0, permutations(i,:))]);
    end
end

function err = Errp(theta,theta0,per)
    % theta = theta.w, theta.g
    % theta.g = theta.g.m , theta.g.s
    k=length(theta.w);
    err = 0;
    for j=1:k
        err = err + theta0.w(j) * H(theta.g{per(j)}.m, theta.g{per(j)}.s, ...
                                    theta0.g{j}.m, theta0.g{j}.s);
    end
    err = err + 0.5 * sum(abs(theta.w(per) - theta0.w));
end

function h = H(m1,s1, m2,s2)
    s = 0.5*(s1 + s2);
    h = (det(s1)^0.25 * det(s2)^0.25) / sqrt(det(s));
    h = h * exp(-1/8 * (m1-m2)'*(s\(m1-m2)));
    h = sqrt(1-h);
end