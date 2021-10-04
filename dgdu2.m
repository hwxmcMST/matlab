function hessian = dgdu2(u,dgdu,model,index,disttype,...
            distpara,para)

n = length(u);
step = max(abs(u./1e4),1e-4);
g1(1:n) = 0;
g2(n,n) = 0;
hessian(n,n) = 0;
order = 0;
x = u2x(u,disttype,distpara);
g = feval(model,x,para,order,index);

for i = 1:n
   tempi = u(i);
   u(i) = u(i) + step(i);
   x = u2x(u,disttype,distpara);
   g1(i) = feval(model,x,para,order,index);
   u(i) = tempi;
   for j = 1:n
      if(i<=j)
         tempi = u(i);
         tempj = u(j);
         u(i) = u(i) + step(i);
         u(j) = u(j) + step(j);
         x = u2x(u,disttype,distpara);
         g2(i,j) = feval(model,x,para,order,index);
         u(i) = tempi;
         u(j) = tempj;
      end
   end
end
for i = 1:n
   for j = 1:n
      if i <= j
         hessian(i,j) = (g2(i,j) - g1(i) - g1(j) + g)/step(i)/step(j);
      else
         hessian(i,j) = hessian(j,i);
      end
   end
end

