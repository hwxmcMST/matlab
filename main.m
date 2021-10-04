% Reliability analysis : SORM-SPA system reliability 
% Hao Wu
% Date: 12/24/2018
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
%----------------------Failure modes with SORM ----------------------------
disp('SORM');
[nr,temp]=size(disttype);
for i=1:nLSF
    [pf_Breitung_t,pf_Tvedt_t,exitFlag,dgdu22_t]=SORM(model(i,:),umpp(i,:),...
        beta(i,:),grad(i,:),disttype,distpara);
    dgdu22{i}=dgdu22_t;
    pf_SORM(i,:)=pf_Breitung_t;
    beta_SORM(i,:)=-norminv(pf_Breitung_t);
    funcall_SORM(i,:)=funcall;
    funcall=0;
    funcall_SORM(i,:)=funcall_FORM(i,:)+funcall_SORM(i,:)-(nr+1);
end

%----------------------Failure modes with S0SPA ---------------------------
disp('SPA SORM');
for i=1:nLSF
    pf_SOSPA_t=spa_sorm(umpp(i,:),grad(i,:),dgdu22{i});
    pf_SOSPA(i,:)=pf_SOSPA_t;
    beta_SOSPA(i,:)= -norminv(pf_SOSPA_t);
    RV_SOSPA(i,:)=1-pf_SOSPA_t;
end
funcall_SOSPA = funcall_SORM;
%----------------------Probability of failure with MCS---------------------
[pfsysMCS,pfMCS]= MCS_SORM(nsimu,nLSF,model,disttype,distpara);

%----------------------System reliability analysis-------------------------
mu_sys_sospa = -beta_SOSPA;
mu_sys_form = -beta;
mu_sys_sorm = -beta_SORM;
sigma = eye(nl,nl);
for i = 1:nl
    for j = i+1:nl
        sigma(i,j) = dot(alpha(i,:),alpha(j,:));
        sigma(j,i) = sigma(i,j);
    end
end
optmvn=optimset('TolFun',1e-6,'MaxFunEvals',1e9); 
rng('default')  % For reproducibility
R_FORM = mvncdf(zeros(nl,1),mu_sys_form,sigma);
R_SORM = mvncdf(zeros(nl,1),mu_sys_sorm,sigma);
R_SOSPA = mvncdf(zeros(nl,1),mu_sys_sospa,sigma);

pf_FORM = 1-R_FORM;
pf_SOSPA = 1-R_SOSPA;
pf_SORM = 1-R_SORM;
error_FORM = abs(pf_FORM-pfsysMCS)/pfsysMCS*100;  
error_SORM = abs(pf_SORM-pfsysMCS)/pfsysMCS*100;
error_SOSPA = abs(pf_SOSPA-pfsysMCS)/pfsysMCS*100;  



disp('-----------------------------------------------------------------------------');
disp('-------------------------------RESULTS---------------------------------------');
disp('-----------------------------------------------------------------------------');
fprintf('pfsysMCS: %.4e\n', pfsysMCS);
fprintf('pf_SOSPA: %.4e\n', pf_SOSPA);
fprintf('pf_SORM: %.4e\n', pf_SORM);
fprintf('pf_FORM: %.4e\n', pf_FORM);
fprintf('error_SOSPA(percentage): %.2e\n', error_SOSPA);
fprintf('error_SORM(percentage): %.2e\n', error_SORM);
fprintf('error_FORM(percentage): %.2e\n', error_FORM);
% fprintf('error_MCS_SEC(percentage): %.2e\n',error_MCS_SEC);
fprintf('funcall_FORM: %d\n', sum(funcall_FORM));
fprintf('funcall_SORM: %d\n', sum(funcall_SORM));
fprintf('funcall_SOSPA: %d\n', sum(funcall_SOSPA));