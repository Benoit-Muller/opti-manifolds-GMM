function question31(d,k,n,scale,incldude_f)
    question31abcdef(d,k,n,scale,incldude_f);
    %question31g(d,k,n,scale)
end

function problem = question31a(d,k,n,scale)
    [mu,sigma,w,xx,n] = makedata(d,k,n,scale,false);
    [u,X,y,Theta] = reparametrize(w,mu,sigma,xx);
    problem.Theta_true = Theta;
    problem.y = y;
    switch 2
        case 1
            problem.M = M_MLE3(d,k);
            problem.cost = @(point) loglikelyhood(point.u,point.X,y);
            problem.egrad = @(point) egrad_l(point.u,point.X,y);
            %problem.egrad = @(point) getApproxGradient(problem, point);
        case 2
            D = diag(randn(d, 1));
            [Q, ~] = qr(randn(d));
            A = Q*D*Q';
            % Create the problem structure.
            problem.M = spherefactory(d);
            % Define the problem cost function and its Euclidean gradient.
            problem.cost  = @(x) -x'*(A*x);
            problem.egrad = @(x) -2*A*x;
    end
end

function  question31abcdef(d,k,n,scale,include_f)
    % a)
    problem = question31a(d,k,n,scale);
    % b,c)
    option = struct("maxiter",Inf, "maxtime",20, "tolgradnorm",1e-5,"verbosity",0);
    endtime = struct("r",[],"c",[]);
    T=30;
    for i=1:T
        option.x0 = problem.M.rand();
        [r.x,r.cost,r.info,r.option] = RGD(problem, option);
        [c.x,c.cost,c.info,c.option] = conjugategradient(problem, option.x0, option);
        endtime.r(i) = r.info(end).time;
        endtime.c(i) = c.info(end).time;
    end
    [r.std, r.mean] = std(endtime.r);
    [c.std, c.mean] = std(endtime.c);
    
    fprintf("Average running times with k=%d, d=%d, n=%d, tolgradnorm=%f :\n", ...
            k,d,n,option.tolgradnorm);
    fprintf("   riemanian GD  : %f +- %f \n",r.mean, r.std)
    fprintf("   congugated GD : %f +- %f \n",c.mean, c.std)

    % d)
    subplot(1,2,1);
    semilogy([r.info.iter], [r.info.gradnorm], '.-');
    hold on;
    semilogy([c.info.iter], [c.info.gradnorm], '.-');
    legend("Riemanian gradient descent", "Conjugated gradient descent")
    xlabel('iteration');
    ylabel('$$\|\nabla l(\Theta_k)\|_{\Theta_k}$$',"Interpreter","latex");
    title("Gradient norm")

    % e)
    subplot(1,2,2);
    plot([r.info.iter], [r.info.cost], '.-');
    hold on;
    plot([c.info.iter], [c.info.cost], '.-');
    legend("Riemanian gradient descent", "Conjugated gradient descent")
    xlabel('iteration');
    ylabel('$$l(\Theta_k)$$',"Interpreter","latex");
    title("Cost")
    sgtitle('Riemann and Conjugated Gradient descent evolution')

    % f)
    if include_f
        [w,mu,sigma,Theta] = deparametrize(c.x.u,c.x.X);
        figure();
        plot(problem.y(1, :), problem.y(2, :), '.', 'MarkerSize', 8)
        hold on
        for j=1:k
            error_ellipse(sigma{j},mu{j});
        end
        title("Data and obtained clusters")
    end
end

function question31g(d,k,n,scale)
    % g)
    option = struct("maxiter",Inf, "maxtime",20, "tolgradnorm",1e-5,"verbosity",0);
    T=30;
    err = zeros(T,1);
    for i=1:T
        problem = question31a(d,k,n,scale);
        option.x0 = problem.M.rand();
        [r.x,r.cost,r.info,r.option] = RGD(problem, option);
        [c.x,c.cost,c.info,c.option] = conjugategradient(problem, option.x0, option);

        [w,mu,sigma,Theta] = deparametrize(x.u,x.X);
        err(i) = Err(Theta, ptob.Theta_true);
    end
    figure();
    h = histogram(err);
    title("Histogram of $$Err(\Theta,\Theta^*)$$ for random data and initial guess")
end








