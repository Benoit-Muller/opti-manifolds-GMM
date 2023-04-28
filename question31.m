function [r,c] = question31(d,k,n,scale,incldude_f)
    [r,c] = question31abcdef(d,k,n,scale,incldude_f);
    [err, h] = question31g(d,k,n,scale);
end

% function problem = question31a(d,k,n,scale)
%     [mu,sigma,w,xx,n] = makedata(d,k,n,scale,false);
%     [u,X,y,Theta] = reparametrize(w,mu,sigma,xx);
%     problem.Theta_true = Theta;
%     problem.y = y;
%     switch 2
%         case 1
%             problem.M = M_MLE3(d,k);
%             problem.cost = @(point) loglikelyhood(point.u,point.X,y);
%             problem.egrad = @(point) egrad_l(point.u,point.X,y);
%             %problem.egrad = @(point) getApproxGradient(problem, point);
%         case 2
%             D = diag(randn(d, 1));
%             [Q, ~] = qr(randn(d));
%             A = Q*D*Q';
%             % Create the problem structure.
%             problem.M = spherefactory(d);
%             % Define the problem cost function and its Euclidean gradient.
%             problem.cost  = @(x) -x'*(A*x);
%             problem.egrad = @(x) -2*A*x;
%     end
% end

function [r,c] = question31abcdef(d,k,n,scale,include_f)
    % a)
    problem = problem_MLE3(d,k,n,scale);
    % b,c)
    T=10;
    option = struct("maxiter",Inf, "maxtime",10, "tolgradnorm",1e-3,"verbosity",0);
    endtime = struct("r",[],"c",[]);
    disp("c)")
    for i=1:T
        fprintf("iter %d/%d\n",i,T)
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
    figure('Position', 10+50*[0 0 14 6])
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
    
    saveas(gcf,sprintf('graphics/q31_evolution_%d_%d_%d.pdf',d,k,n))

    % f)
    if include_f
        [w,mu,sigma,Theta] = deparametrize(c.x.u,c.x.X);
        figure();
        plot(problem.y(1, :), problem.y(2, :), '.', 'MarkerSize', 8)
        hold on
        for j=1:k
            if all(eig(sigma{j})>0)
                error_ellipse(sigma{j},mu{j});
            else
                warning("sigma non positive definite")
            end
        end
        title("Data and obtained clusters")
        saveas(gcf,sprintf('graphics/q31_klusters_%d_%d_%d.pdf',d,k,n))
    end
end

function [err, h] = question31g(d,k,n,scale)
    % g)
    option = struct("maxiter",Inf, "maxtime",20, "tolgradnorm",1e-3,"verbosity",0);
    T=20;
    err = zeros(T,1);
    disp("g)")
    for i=1:T
        fprintf("iter %d/%d\n",i,T)
        problem = problem_MLE3(d,k,n,scale);
        option.x0 = problem.M.rand();
        [r.x,r.cost,r.info,r.option] = RGD(problem, option);
        [c.x,c.cost,c.info,c.option] = conjugategradient(problem, option.x0, option);

        err(i) = Err(c.x, problem.theta);
    end
    figure();
    h = histogram(err);
    title("Histogram of $$Err(\Theta,\Theta^*)$$ for random data and initial guess")
    saveas(gcf,sprintf('graphics/q31_hist_%d_%d_%d.pdf',d,k,n))
end








