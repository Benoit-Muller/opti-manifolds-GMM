function g = egrad_l(u,X,y)
    k=length(u);
    [~,n]=size(y);
    Q = zeros(n,k);
    for i=1:n
        Q(i,:) = cellfun(@(X)q(X, y(:,i)), X);
    end
    p = Q * u.^2;
    g.S = -2*u.*(Q'*(1./p));
    g.P=cell(k,1);
    for j=1:k
        g.P{j} = 0;
        for i=1:n
            g.P{j} = g.P{j} + Q(i,j) / p(i) * (X{j} - y(:,i)*y(:,i)');
        end
         g.P{j} = g.P{j} * u(j)^2 / 2;
    end
end
