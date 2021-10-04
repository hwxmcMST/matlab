function [pf,umpp,xmpp,beta,gradient,gmpp] = FORM(model,disttype,...
    distpara,para,derivative,index)

[n,temp]=size(disttype);
u0(1,n) = 0;
[umpp,xmpp,beta,gradient,gmpp] = MPP_HLRF(u0,model,...
 disttype,distpara,para,derivative,index);
pf = normcdf(-beta);
