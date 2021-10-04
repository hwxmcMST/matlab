function [pfs,pfs_low,pfs_high,pf] = MCS(nsimu,model,disttype,distpara,...
    para,nc)
% MCS for ssytem reliability analysis
% input:
% nsimu: number of the simulation point;
% model includes limit state function and related parameters;
% distype is the distribution type of the random variables in this model;
% distpara is the distribution parameter of the random variables in this;
% para is the deterministic parameter in this model;
% nc
% output:
% pfs is the probability of failure for system;
% pfs_low is the low boundary of confidence interval;
% pfs_high is the high boundary of confidence interval;
% pf is the probability of failure of component in the system
% pf = Pr( g < 0)
[nr,col] = size(distpara);
rng('default');
rng(0);
iteration = 1;
n = nsimu;
if nsimu > 1e6
    n=1e6;
    iteration = nsimu/1e6;
end
if nsimu > 1e6 && mod(nsimu,1e6)~=0
    display('Wrong number of simulation');
    return;
end

x = zeros(n,nr);

for k = 1:iteration
    for i = 1:nr
        x(1:n,i)=sample(n,disttype(i,:),distpara(i,1:col));
    end
    for j = 1:nc 
        g(:,j) = feval(model(j,:),x);
    end
    nf(k,:) = sum(g < 0);
    if nc ==1
        nfs = nf;
    end
    if nc > 1
        indicator = g < 0; % 0 = failure
        nfs(k) = sum(sum(indicator')>0);
    end
end
pf = sum(nf)/nsimu;
pfs = sum(nfs)/nsimu;

%Confidence interval
alpha = 1-0.95;
z = norminv(1-alpha/2);
pfs_low = pfs*(1-z*((1-pfs)/(nsimu*pfs))^0.5);
pfs_high = pfs*(1+z*((1-pfs)/(nsimu*pfs))^0.5);

function x=sample(n,type,p)
if strcmp(deblank(type),'Normal') | strcmp(deblank(type),'norm')
    x=normrnd(p(1),p(2),n,1);
end
if strcmp(deblank(type),'Lognormal') | strcmp(deblank(type),'logn')
    b=(log((p(2)/p(1))^2+1))^0.5;
    a=log(p(1))-0.5*b^2;
    x=lognrnd(a,b,n,1);
      
%       b = p(2);
%       a = p(1);
%       x = lognrnd(a,b,n,1);
end
if strcmp(deblank(type),'Weibull') | strcmp(deblank(type),'weib')
    x=weibrnd(p(1),p(2),n,1);
end
if strcmp(deblank(type),'Exponential') | strcmp(deblank(type),'exp')
    x=exprnd(p(1),n,1);
end
if strcmp(deblank(type),'Uniform') | strcmp(deblank(type),'unif')
    x=unifrnd(p(1),p(2),n,1);
end

if strcmp(deblank(type),'Noncenchi') | strcmp(deblank(type),'ncx2')
    x=ncx2rnd(p(1),p(2),n,1);
end


if strcmp(deblank(type),'Extreme1') | strcmp(deblank(type),'ext1')
    z=unifrnd(0,1,n,1);
    alfa=1/6^0.5*pi/p(2);
    mu=p(1)+psi(1)/alfa;
    x=-log(-log(z))/alfa+mu;
end
if strcmp(deblank(type),'tri')
    z=unifrnd(0,1,n,1);
    a=p(1);b=p(2);
    x=(b-a)*z.^0.5+a;
end


