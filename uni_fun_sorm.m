function [k,kp,kpp,c,sigma]=uni_fun_sorm(umpp,grad,dgdu22,t)
%Retunr CGF of the function
umpp=umpp';        % change to column vector
grad=grad';

nr=length(umpp);

a_original=0.5*dgdu22;
b_original=grad-dgdu22*umpp;
c_original=0.5*umpp'*dgdu22*umpp-grad'*umpp;

[P,A_bar]=eig(a_original);
b_bar=P'*b_original;

ki=[];kip=[];kipp=[];
for i=1:nr
    a(i)=c_original/nr;
    b(i)=b_bar(i);
    c(i)=A_bar(i,i);
    
    if c(i)>0
        e(i)=a(i)-b(i)^2/4/c(i);
        mu(i)=0.5*b(i)/c(i)^0.5;
        sigma(i)=c(i)^0.5;
        lamda=(mu(i)/sigma(i))^2;
        ki(i)=e(i)*t+lamda*sigma(i)^2*t/(1-2*sigma(i)^2*t)...
            -0.5*log(1-2*sigma(i)^2*t);
        kip(i)=e(i)+lamda*sigma(i)^2/(1-2*sigma(i)^2*t)+...
            2*lamda*sigma(i)^4*t/(1-2*sigma(i)^2*t)^2+sigma(i)^2/(1-2*sigma(i)^2*t);
        kipp(i)=2*sigma(i)^4*(-1+2*sigma(i)^2*t-2*lamda)/(-1+2*sigma(i)^2*t)^3;
    end
    
    if c(i)<0
        e(i)=a(i)-b(i)^2/4/c(i);
        mu(i)=-0.5*b(i)/(-c(i))^0.5;
        sigma(i)=(-c(i))^0.5;
        lamda=(mu(i)/sigma(i))^2;
        ki(i)=e(i)*t-lamda*sigma(i)^2*t/(1+2*sigma(i)^2*t)...
            -0.5*log(1+2*sigma(i)^2*t);
        kip(i)=e(i)-lamda*sigma(i)^2/(1+2*sigma(i)^2*t)+...
            2*lamda*sigma(i)^4*t/(1+2*sigma(i)^2*t)^2-sigma(i)^2/(1+2*sigma(i)^2*t);
        kipp(i)=2*sigma(i)^4*(2*lamda+1+2*sigma(i)^2*t)/(1+2*sigma(i)^2*t)^3;
    end
    if abs(c(i))<1e-9
        sigma(i)=0;
        ki(i)=a(i)*t+0.5*b(i)^2*t^2;
        kip(i)=a(i)+b(i)^2*t;
        kipp(i)=b(i)^2;
    end

end

k = sum(ki);
kp= sum(kip);
kpp = sum(kipp);

    
    