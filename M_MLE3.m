function M = M_MLE3(d,k)
    % manifold S^(k-1) * (P_(d,br))^k

    Sphere = spherefactory(k);

    spd_br = spd_br_factory(d);
    prod_spd_br = powermanifold(spd_br, k);
    prod_spd_br = rmfield(prod_spd_br,{'exp'});

    M = productmanifold(struct('S',Sphere, 'P',prod_spd_br));