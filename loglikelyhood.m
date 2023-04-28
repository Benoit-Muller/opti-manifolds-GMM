function l = loglikelyhood(u,X,y)
% loglikelyhood l
    [~,n] = size(y);
    l = 0;
    for i=1:n
        l = l - log(cellfun(@(Y)q(Y,y(:,i)),X)' * u.^2);
    end
end

% function X = proj_symbr_matrix(X)
%     X = X - X(end,end) .* X(end,:) .* X(:,end);
% end
% proj_symbr = @(X) cellfun(proj_symbr_matrix,X)

% function [gu,gX] = rgrad_l(u,X,y)
%     [gu,gX] = egrad_l(u,X,y)
%     gX = proj_symbr(gX)
% end

