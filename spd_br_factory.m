function M = spd_br_factory(d)
    % P_(br,d+1): spd matrices X with X_(d+1,d+1)=1
    symm = @(X) .5*(X+X');  
    M.name = @() sprintf(['Symmetric positive definite geometry of' ...
        ' %dx%d matrices with unitary bottom rigth entry'], d+1, d+1);
    M.dim = @() d*(d+3)/2;
    % Helpers to avoid computing full matrices simply to extract their trace
    vec  = @(A) A(:);
    trAB = @(A, B) vec(A')'*vec(B);  % = trace(A*B)
    trAA = @(A) sqrt(trAB(A, A));    % = sqrt(trace(A^2))

    function A = normalize_br(A)
        if abs(A(end,end)) < 1e-100
            warning("division by zero in br normalization");
        end
        A = A/A(end,end);
    end
    
    M.inner = @(X, eta, zeta) trAB(X\eta, X\zeta); %Question 17.a)
    M.norm = @(X, eta) real(trAA(X\eta));
    
    function Y = retraction_spd(X, eta, t)
        eta = symm(eta);
        if nargin < 3
         teta = eta;
        else
            teta = t*eta;
        end
        Y = symm(X + teta + .5*teta*(X\teta));
        % The symm() call is mathematically unnecessary but numerically
    end
    
    M.retr = @retraction;
    function Y = retraction(X, eta, t)
        % Question 7.b)
        Y = normalize_br(retraction_spd(X, eta, t));
    end
    
    M.proj = @projection;
    function eta = projection(X, eta)
        eta = symm(eta);
        eta = eta - (eta(end,end)*X(:,end))*X(end,:);
    end

    M.tangent = M.proj;

    M.egrad2rgrad = @egrad2rgrad;
    function eta = egrad2rgrad(X, eta)
        % Question 7.c), formula of Question 16
        eta = M.proj(X,X*symm(eta)*X);
    end

    M.ehess2rhess = @ehess2rhess;
    function Hess = ehess2rhess(X, egrad, ehess, eta)
        Hess = M.proj(X, X*symm(ehess)*X + eta*symm(egrad)*X);
    end
    
    M.rand = @random;
    function X = random()
        D = diag(1+rand(d, 1));
        [Q, R] = qr(randn(d)); %#ok
        sigma = Q*D*Q';
        mu = randn(d,1);
        X = [sigma + mu*mu' , mu;
             mu'        , 1];
    end

    M.randvec = @randomvec;
    function eta = randomvec(X)
        eta = symm(randn(d+1));
        eta(end,end) = 0;
        eta = eta / M.norm(X, eta);
    end

    M.lincomb = @matrixlincomb;
    M.transp = @(X1, X2, eta) eta;
    M.zerovec = @(X) zeros(d+1);
    M.vec = @(X, U) U(:);
    M.mat = @(X, u) reshape(u, d+1, d+1);
    M.vecmatareisometries = @() false;

end