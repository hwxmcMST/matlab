function dgdu22=hessian(model,umpp,disttype,distpara)
%2nd-order derivatives (finite difference)
%Xiaoping Du (10/16/09)
global funcall;
x=u2x(umpp,disttype,distpara);
gmpp=feval(model,x);

% step1(1:length(umpp))=1.e-6;
step1=min(8e-7,max(abs(umpp/100),1.0e-5));
step2=step1;
[n,m]=size(disttype);
gDel2=zeros(n,n);
for i=1:n
   u=umpp;u(i)=u(i)+step1(i);x=u2x(u,disttype,distpara);
   gDel1(i)=feval(model,x);
   for j=1:n
      if(i<=j)
         u=umpp;u(i)=u(i)+step1(i);u(j)=u(j)+step2(j);
         x=u2x(u,disttype,distpara);
         gDel2(i,j)=feval(model,x);
      end
   end
end
dgdumpp=(gDel1-gmpp)./step1;
for i=1:n
   for j=1:n
      if i<=j
         dgdu22(i,j)=(gDel2(i,j)-gDel1(i)-gDel1(j)+gmpp)/step1(i)/step2(j);
      else
         dgdu22(i,j)=dgdu22(j,i);
      end
   end
end
