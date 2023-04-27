function [x, cost, info, options] = RGD(problem, options)
% options.x0      : initial guess (optional)
%        .maxiter : (optional)
%        .tol     : (optional)
%        .tau     : (optional)
%        .r       : (optional)
    if ~ isfield(options,"x0")
        x0 = problem.M.rand();
    end
    if ~ isfield(options,"maxiter")
        options.maxiter = 100;
    end
    if ~ isfield(options,"tol")
        options.tol = 1e-5;
    end
    if ~ isfield(options,"tau")
        options.tau = 1/2;
    end
    if ~ isfield(options,"r")
        options.r = 1e-4;
    end

    x = x0;
    iter = 0;
    while iter < options.maxiter
        [cost,g] = getCostGrad(problem,x);
        if problem.M.norm(x,g) < options.tol
            break
        end
        line = @(alpha) problem.cost(M.retr(x,alpha*g));
        alpha = backtracking(problem.M,cost,g,line,options.tau,options.r);
        s = -alpha * g;
        x = M.retr(x,s);
        iter = iter + 1;
    end
    info.iter = iter;
    info.gradnorm = problem.M.norm(x,g);
end

function alpha = backtracking(M,x,f,g,line,tau,r)
    alpha = 1;
    maxiter = 100;
    iter=1;
    while f - line(alpha) < r * M.norm(x,g)^2 && iter <= maxiter
        alpha = tau * alpha;
        iter = iter + 1;
    end
end