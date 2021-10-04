function [umpp,xmpp,beta,gradient,gmpp]=MPP_HLRF(u0,model,disttype,distpara,...
    para,derivative,index)

% Initialization
funcall_FORM = 0;
[n,temp] = size(disttype);
convergence = 0;
tol_beta=1e-4; iter_max=300; iter = 0; funcall = 0;
u = u0; % starting point
beta = norm(u);
order = 1;
while ~convergence
    iter = iter + 1;
    [g,gradient] = gatu(model,u,disttype,distpara,para,order,index,derivative);
    beta0 = beta;
    beta = beta + g / norm(gradient);
    alpha= gradient / norm(gradient);
    u = -beta * alpha;
    % Check convergence
    if abs(beta - beta0) < tol_beta 
        convergence =1;
    end
    if iter >= iter_max
        disp('Maximum number of iterations is reached for MPP search.')
        return
    end
%     funcall_FORM = funcall_FORM + (n + 1);
end
umpp = u;
xmpp = u2x(u,disttype,distpara);
gmpp =g;
