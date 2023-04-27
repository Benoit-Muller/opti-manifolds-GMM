function [x, cost, info, option] = RGD(problem, option)
% option.x0      : initial guess (optional)
%        .maxiter : (optional)
%        .tol     : (optional)
%        .tau     : (optional)
%        .r       : (optional)
    tic
    if ~ isfield(option,"x0")
        option.x0 = problem.M.rand();
    end
    if ~ isfield(option,"maxiter")
        option.maxiter = inf;
    end
    if ~ isfield(option,"tolgradnorm")
        option.tolgradnorm = 1e-5;
    end
    if ~ isfield(option,"tau")
        option.tau = 1/2;
    end
    if ~ isfield(option,"r")
        option.r = 1e-4;
    end
    if ~ isfield(option,"maxtime")
        option.maxtime = 10;
    end

    info=struct('iter',[],'gradnorm',[],'time',[],'cost',[],'alpha',[]);
    x = option.x0;
    iter = 0;
    while iter < option.maxiter && toc < option.maxtime
        [cost,g] = getCostGrad(problem,x);
        gradnorm = problem.M.norm(x,g);
        if gradnorm < option.tolgradnorm
            break
        end
        line = @(alpha) problem.cost(problem.M.retr(x,g,-alpha));
        alpha = backtracking(problem.M,x,cost,g,line,option.tau,option.r);
        info(iter+1) = struct('iter',iter,'gradnorm',gradnorm, ...
                              'time',toc,'cost',cost,'alpha',alpha); 
        x = problem.M.retr(x, g,-alpha);
        iter = iter + 1;
    end
end

function alpha = backtracking(M,x,f,g,line,tau,r)
    gnorm = M.norm(x,g);
    alpha = 1/tau;
    maxiter = 100;
    iter=1;
    while f - line(alpha) < r * alpha * gnorm^2 && iter <= maxiter
        alpha = tau * alpha;
        iter = iter + 1;
    end
end