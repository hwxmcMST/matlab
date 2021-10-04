function [umpp,xmpp,beta,grad,gmppx,flag]=MPP_SQP(u0,model,disttype,distpara)
%%Use SQP for MPP search
%10/17/2009
%Failure region: g-gmpp<0;
n=length(u0);
UB(1:n)=15;
LB=-UB;
% option=optimset('display','off','MaxFunEvals',1000,'Tolfun',1e-6,'Algorithm','active-set');
option=optimset('display','off','MaxFunEvals',1000);
[umpp,~,flag]=fmincon('mpp_obj',u0,[],[],[],[],LB,UB,'mpp_con',option,...
    model,disttype,distpara);
beta=norm(umpp);
xmpp=u2x(umpp,disttype,distpara);
g=feval(model,xmpp);
gmppx=g;
stepsize=1e-3;
for i=1:n
        u_temp=umpp;u_temp(i)=u_temp(i)+stepsize;
        x_temp=u2x(u_temp,disttype,distpara);
        g_temp=feval(model,x_temp);g_temp=g_temp;
        grad(i)=(g_temp-g)/stepsize;
end
