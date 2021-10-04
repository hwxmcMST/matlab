function [g,gradient,hessian] = gatu(model,u,disttype,distpara,...
    para,order,index,derivative)

n = length(u);
hessian(n,n) = 0;
x = u2x(u,disttype,distpara);
if order == 0
    g = feval(model,x,para,index); 
end
if order==1 || order == 2
    g = feval(model,x,para,index); 
    if order == 1 && derivative == 0
        gradient = dgdu(u,g,model,index,disttype,distpara,...
            para);
    end
    if order == 2 && derivative == 0
        gradient = dgdu(u,g,model,index,disttype,distpara,...
            para);
        hessian = dgdu2(u,gradient,model,index,disttype,...
            distpara,para);
    end
    if order == 1 && derivative == 1
        dXdU = dxdu(u,disttype,distpara);
        gradient = dgdx .* dXdU;
    end
    if order == 2 && derivative == 1
        dXdU = dxdu(u,disttype,distpara);
        gradient = dgdx .* dXdU;
        for i = 1:n
            for j = 1:n
                hessian(i,j) = dgdx2(i,j) * dXdU(i) * dXdU(j);
            end
        end
    end
end


