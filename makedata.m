function [mu,Sigma,w,X] = makedata(d,k,n,scale,display)
    % create parameters w, mu_j, Sigma_j for the ground truth GMM
    for j = 1 : k
        [mu{j}, Sigma{j}] = randomgaussian(d, scale);
    end
    w = randfixedsum(k,1,1,0,1); % select w from simplex
    
    % sample and plot data points x_1, ..., x_n from the ground truth GMM
    % aggregate the samples into the matrix X; the columns of X are the samples x_i = X(:, i)
    X = samplegaussian(mu{1}, Sigma{1}, round(n*w(1)));
    if d==2 && display
        plot2Ddata(X);
        hold on;
    end
    for j = 2 : k
        tempX = samplegaussian(mu{j}, Sigma{j}, round(n*w(j)));
        if d==2 && display
            plot2Ddata(tempX);
        end
        X = [X, tempX]; % aggregate samples into matrix X
    end
    n = size(X, 2); % true number of samples 
end
 
%%%  Helpers:  %%%

function plot2Ddata(Y)
    plot(Y(1, :), Y(2, :), '.', 'MarkerSize', 8);
end

% generate random vector mu and positive definite matrix Sigma
function [mu, Sigma] = randomgaussian(d, scale)
    mu = randn(d, 1);
    Y = randn(d, d);
    Sigma = scale^2*Y.'*Y; 
    % here, you can also add the term " + const*trace(scale^2*Y.'*Y)*eye(d) "
    % to reduce the eccentricity of the Gaussian
    Sigma = Sigma + trace(scale^2*Y.'*Y)*eye(d);
end

% sample n points from the Gaussian with center mu and covariance Sigma
function Y = samplegaussian(mu, Sigma, n)
    d = numel(mu);
    R = chol(Sigma);
    Y = R'*randn(d, n) + mu;
end
