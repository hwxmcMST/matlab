function dXdU = dxdu(u,disttype,distpara)
% Find derivatives of X wrt U
% 11/12/2017
n = length(u);
dXdU(1:n) = 0;
for i = 1:n
    if strcmp(deblank(disttype(i,:)),'Normal') || strcmp(deblank(disttype(i,:)),'norm')
        dXdU(i) = distpara(i,2);
    end
    if strcmp(deblank(disttype(i,:)),'Lognormal') || strcmp(deblank(disttype(i,:)),'logn')
       % distpara(i,1)-mean of x and distpara(i,2)-std of x
       % a - mean of log(x); b - variance of log(x)
       b = (log((distpara(i,2)/distpara(i,1))^2+1))^0.5;
       a = log(distpara(i,1))-0.5*b^2;
       dXdU(i) = normpdf(u(i)) / lognpdf(logninv(normcdf(u(i)),a,b),a,b);
   end
   if strcmp(deblank(disttype(i,:)),'Weibull') || strcmp(deblank(disttype(i,:)),'weib')
       dXdU(i) = normpdf(u(i)) / wblpdf(wblinv(normcdf(u(i)),distpara(i,1),distpara(i,2)),...
           distpara(i,1),distpara(i,2));
   end
   if strcmp(deblank(disttype(i,:)),'Exponential') || strcmp(deblank(disttype(i,:)),'exp')
       dXdU(i) = normpdf(u(i)) / exppdf(expinv(normcdf(u(i)),distpara(i,1)),...
           distpara(i,1));
   end
   if strcmp(deblank(disttype(i,:)),'Extreme1') | strcmp(deblank(disttype(i,:)),'ext1')
   %Maximum, f(x)=1/b*exp(-(x-a)/b)*exp(-exp(-(x-a)/b));   
        p1 = distpara(i,1);
        p2 = distpara(i,2);
        z = normcdf(u(i));
        x = evmaxinv(z,p1,p2);
        pdf_x = evmaxpdf(x,p1,p2);
        dXdU(i) = normpdf(u(i)) / pdf_x;
    end
end
