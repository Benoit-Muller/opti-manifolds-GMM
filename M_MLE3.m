function problem = problem_MLE3(d,k,n,scale)
    % manifold S^(k-1) * (P_(d,br))^k

    Sphere = spherefactory(k);

    spd_br = spd_br_factory(d);
    prod_spd_br = powermanifold(spd_br, k);
    prod_spd_br = rmfield(prod_spd_br,{'exp'});

   problem.M = productmanifold(struct('u',Sphere, 'X',prod_spd_br));

    display=false;
    [mu,sigma,w,xx,problem.n] = makedata(d,k,n,scale,display);
    [problem.u,problem.X,problem.y,problem.Theta] = reparametrize(w,mu,sigma,xx);

    problem.cost = @(point) loglikelyhood(point.u,point.X,y);
    problem.egrad = @(point) egrad_l(point.u,point.X,y);
end