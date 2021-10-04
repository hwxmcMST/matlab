function x=u2x(u,type,p)
% Transform a standard normal variable to a non-normal variable
%Xiaoping Du (09/20/09)
n=length(u);k=0;
for i=1:n
    if strcmp(deblank(type(i,:)),'Normal') | strcmp(deblank(type(i,:)),'norm')
        x(i)=p(i,1)+u(i)*p(i,2);
    end
    if strcmp(deblank(type(i,:)),'Lognormal') | strcmp(deblank(type(i,:)),'logn')
       %p(1)-mean of x and p(2)-std of x  (p1 & p2 are lognormal parameter)
       %a-mean of log(x); b-variance of log(x) (a & b are normal parameter)
       b=(log((p(i,2)/p(i,1))^2+1))^0.5;
       a=log(p(i,1))-0.5*b^2;
       x(i)=logninv(normcdf(u(i)),a,b);
   end
   if strcmp(deblank(type(i,:)),'Weibull') | strcmp(deblank(type(i,:)),'weib')
       mu = p(i,1);
       sigma = p(i,2);
       f=@(k)sigma^2/mu^2-gamma(1+2/k)/(gamma(1+1/k)^2)+1;
       k0 = (sigma/mu).^-1.086;
       k = fzero(f,k0);
       lam=mu/gamma(1+1/k); 
       x(i)=wblinv(normcdf(u(i)),lam,k);
   end
   if strcmp(deblank(type(i,:)),'Exponential') | strcmp(deblank(type(i,:)),'exp')
       x(i)=expinv(normcdf(u(i)),p(i,1));
   end
   if strcmp(deblank(type(i,:)),'Uniform') | strcmp(deblank(type(i,:)),'unif')
       x(i)=unifinv(normcdf(u(i)),p(i,1),p(i,2));
   end
   if strcmp(deblank(type(i,:)),'Chisquare') | strcmp(deblank(type(i,:)),'chi2')
       x(i)=chi2inv(normcdf(u(i)),p(i,1));
   end
   if strcmp(deblank(type(i,:)),'Extreme1') | strcmp(deblank(type(i,:)),'ext1')
%Maximum, f(x)=1/b*exp(-(x-a)/b)*exp(-exp(-(x-a)/b));       
       z=normcdf(u(i));
       b=6^0.5*p(i,2)/pi;
       a=p(i,1)+Psi(1)*b; 
       z=max(z,10e-60);
       x(i)=a-log(-log(z))*b;
   end
   if strcmp(deblank(type(i,:)),'tri') 
       a=p(i,1);b=p(i,2);
       x(i)=(b-a)*normcdf(u(i))^0.5+a;
   end
   if strcmp(deblank(type(i,:)),'clearance') & k~=i
       r=p(i,1);
       if r<1.0e-20
           x(i)=0;x(i+1)=0;
       end
       if r>=1.0e-20
          [x(i),f]=fminbnd('myfun_obj',-r,r,[],u(i),r);
%          [x(i+1),f]=fminbnd('myfun_obj',-r,r,[],u(i+1)   ,r);
           if f>1e-6
              display('failed');
             pause
           end
           x(i+1)=(2*normcdf(u(i+1))-1)*(r^2-x(i)^2)^0.5;
%           x(i)=(2*normcdf(u(i))-1)*(r^2-x(i+1)^2)^0.5;
       end
       k=i+1;
   end
end
