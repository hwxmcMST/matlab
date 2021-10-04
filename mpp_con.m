function [con,ceq]=mpp_con(u,model,disttype,distpara)
%Constraint function for MPP search
x=u2x(u,disttype,distpara);
g=feval(model,x);
con=[];
ceq=g;
