function g = egrad_l(u,X,y)
    k=length(u);
    [~,n]=size(y);
    Q = zeros(n,k);
    for i=1:n
        Q(i,:) = cellfun(@(X)q(X, y(:,i)), X);
    end
    p = Q * u.^2;
    g.u = -2*u.*(Q'*(1./p));
    g.X=cell(k,1);
    for j=1:k
        g.X{j} = 0;
        for i=1:n
            g.X{j} = g.X{j} + Q(i,j) / p(i) * (X{j} - y(:,i)*y(:,i)');
        end
         g.X{j} = g.X{j} * u(j)^2 / 2;
    end
end
