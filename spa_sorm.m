function pf = spa_sorm(umpp,grad,dgdu22)
%Return pf
%Xiaoping Du, 03/14/10

%Solve saddlepoint
t0=0;
options=optimset('Algorithm','interior-point','display','iter','MaxFunEvals',10000);
t_sorm = fmincon('solve_sp_obj_sorm',t0,[],[],[],[],[],[],'solve_sp_con_sorm',options,...
    umpp,grad,dgdu22);

[k,kp,kpp,~,~]=uni_fun_sorm(umpp,grad,dgdu22,t_sorm);

gmpp=0;
w=sign(t_sorm)*(2*(t_sorm*gmpp-k))^0.5;
v=t_sorm*(kpp)^0.5;
pf=normcdf(w)+normpdf(w)*(1/w-1/v);

