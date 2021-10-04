function [pf_Breitung,pf_Tvedt,exitFlag,dgdu22]=SORM(model,umpp,...
    beta,dgdu1,disttype,distpara)
%Second Order Reliability Method (SORM)
%Xiaoping Du (10/16/2009)
global funcall;
[n,m]=size(disttype);
%exitFlag: 1 - successful, -3 - sigularity 
exitFlag=1;
%1st-order and 2nd-order derivatives
dgdu22=hessian(model,umpp,disttype,distpara);

%Gram-Schmidt Othogonalization
r0=eye(n,n);
r0(n,:)=-dgdu1/norm(dgdu1);
r=zeros(n,n);
r(n,:)=r0(n,:);
for k=n-1:-1:1
    r(k,:)=r0(k,:);
    for j=n:-1:k+1
        r(k,:)=r(k,:)-(dot(r(j,:),r0(k,:))/norm(r(j,:))^2)*r(j,:);
    end
    r(k,:)=r(k,:)/norm(r(k,:));
end
a=r*dgdu22*r'/norm(dgdu1);
A=a(1:n-1,1:n-1);
[v,d]=eig(A);

%Calculate probability
B1=1; B2=1; B3=1;
for i=1:n-1
    %Check sigularity
    if (1+beta*d(i,i))<0
        exitFlag=-3;
        break;
    end
    B1=((1+beta*d(i,i))^(-0.5))*B1;
    B2=((1+(beta+1)*d(i,i))^(-0.5))*B2;
    B3=((1+(beta+1i)*d(i,i))^(-0.5))*B3;
end
A1=normcdf(-beta)*B1;
A2=(beta*normcdf(-beta)-normpdf(beta))*(B1-B2);
A3=(beta+1)*(beta*normcdf(-beta)-normpdf(beta))*(B1-real(B3));
pf_Breitung=A1;
pf_Tvedt=A1+A2+A3;