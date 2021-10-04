function g1=exp1_1(x,~,para,index)
global funcall;
funcall=funcall+1;
x1=x(:,1);
x2=x(:,2);

% g < 0 is a failure
g1= (x1.*x2-5);