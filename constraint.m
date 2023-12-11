function [c,ceq]=constraint(x)
in=length(x)-7;
in=in/2;
%control points
P=zeros(2,in+2);
P(1:2,1)=(0);%���
P_end = [0.03;0.015];
P(1:2,in+2)=P_end;%�յ�

for i=1:in
    P(1,i+1)=x(2*i-1);%�м���Ƶ�x����
    P(2,i+1)=x(2*i);%�м���Ƶ�y����
end

nn=0;
c=zeros(1,in);
for i=1:in
    nn=nn+1;
    L1=(P(1,i)-P(1,i+1))^2+(P(2,i)-P(2,i+1))^2;
    L2=(P(1,i+2)-P(1,i+1))^2+(P(2,i+2)-P(2,i+1))^2;
    L3=(P(1,i)-P(1,i+2))^2+(P(2,i)-P(2,i+2))^2;
    l1=sqrt(L1);
    l2=sqrt(L2);
    pos=(L1+L2-L3)/(2*l1*l2);
    angle=acos(pos);
    c(nn)=real(pi/3-angle);
end
ceq=[];
end
