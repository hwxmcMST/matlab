function gradient = dgdu(u,g,model,index,...
    disttype,distpara,para)

% Xiaoping Du, 01/31/2018
n = length(u);
step = max(abs(u./1e3),1e-4);
gradient(1:n) = 0;
order = 0;
for i = 1:n
    temp = u(i);
    u(i) = u(i) + step(i);
    x = u2x(u,disttype,distpara);
    gNew = feval(model,x,para,index);
    gradient(i) = (gNew - g) / step(i);
    u(i) = temp;
end

