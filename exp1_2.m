function g2=exp1_2(x, ~,para,index)
global funcall;
funcall=funcall+1;
x1=x(:,1);
x2=x(:,2);

g2 = x1.^2+x2.^2+8*x1-16*x2+40;