clear; close all; warning off; clc;
exID=5;  

[model,disttype,distpara,nLSF]=exp_in(exID);
para=0;
index=1;
derivative=0;
nsimu=1e7;
gmpp=0;
global funcall;
funcall=0;
%----------------------Failure modes with FORM ----------------------------
disp('FORM');
nl=nLSF;
for i=1:nLSF
    [pf_FORM_t,umpp_t,xmpp_t,beta_t,grad_t,gmppx_t] = FORM(model(i,:),disttype,...
    distpara,para,derivative,index);
    pf_FORM(i,:)=pf_FORM_t;
    grad(i,:)=grad_t;
    alpha(i,:)=grad_t/norm(grad_t);
    beta(i,:)=beta_t;
    umpp(i,:)=umpp_t;
    xmpp(i,:)=xmpp_t;
    gmppx(i,:)=gmppx_t;
    funcall_FORM(i,:)=funcall;
    funcall=0;
end