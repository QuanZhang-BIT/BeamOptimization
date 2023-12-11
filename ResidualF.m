function [Delta_F]=ResidualF(x)
tic;
Target_F = load('E:\QuanZhang\0-0_3.0-1.5cm(disp0.8(0.5-1.0))thickness1mm\Target_F.txt');
in=length(x)-7;
in=in/2;

model=mphload('model_test.mph');

%control points
P=zeros(2,in+2);
P(1:2,1)=(0);%起点
P_end = [0.03;0.015];
P(1:2,in+2)=P_end;%终点

for i=1:in
    P(1,i+1)=x(2*i-1);%中间控制点x坐标
    P(2,i+1)=x(2*i);%中间控制点y坐标
end
kk=0;
for i=1:in-1
    for j=i+2:in+1
        kk=kk+1;
        r0=(P(:,i)+P(:,i+1))/2-P(:,j);
        rs=P(:,i)-P(:,i+1);
        R90=[0 -1;1 0];
        d0=R90*rs/norm(rs);
        n00=dot(r0,d0)*d0/norm(dot(r0,d0)*d0);
        r1=P(:,i)-P(:,j);
        r2=P(:,i+1)-P(:,j);
        d1=R90*r1/norm(r1);
        d2=R90*r2/norm(r2);
        n1=-sign(dot(d1,rs))*d1;
        n2=sign(dot(d2,rs))*d2;
        mu=det([r0 rs]);
        theta=10/180*pi;
        n10=RotateM(sign(mu)*theta)*n1;
        n20=RotateM(-sign(mu)*theta)*n2;
        f0=dot(n00,P(:,j+1)-P(:,i));
        f1=dot(n10,P(:,j+1)-P(:,i));
        f2=dot(n20,P(:,j+1)-P(:,i+1));
        c(kk)=min([f0 f1 f2]);
    end
end

if max(c)>0
    Delta_F=100;
else
n = 3; %curve's degree n = p + 1 where p is the degree of the curve
%knot vector
t = [0 0 0 0.3 0.5 0.7 1 1 1];
%weight vector
%w = [1 1 1 1 1];
w = x(end-5:end);
%w = [1 x(end-4:end)];
%call of the function
C = nurbsfun(n, t, w, P);
nOrder=size(C,2);
% dis=zeros(nOrder-1,1);
% mini_dis=0;
for ii=1:nOrder
    strX_value=C(1,ii);strY_value=C(2,ii);
    model.geom('geom1').feature('pol1').setIndex('table', strX_value, ii-1, 0);
    model.geom('geom1').feature('pol1').setIndex('table', strY_value, ii-1, 1);
end
model.geom('geom1').run('pol1');
model.geom('geom1').run;
bb=x(end-6)/10000;
model.param.set('b',bb);
model.param.set('ht',0.01);

coord = P_end';
NodeNum = mphselectcoords(model,'geom1',coord','point');
model.component('comp1').cpl.create('aveop1', 'Average');
model.component('comp1').cpl('aveop1').set('axisym', true);
model.component('comp1').cpl('aveop1').selection.geom('geom1', 0);
model.component('comp1').cpl('aveop1').selection.set([NodeNum]);

model.component('comp1').variable.create('var1');
model.component('comp1').variable('var1').set('Fx', 'aveop1(beam.RFx)');
model.component('comp1').variable('var1').set('Fy', 'aveop1(beam.RFy)');
model.component('comp1').variable('var1').set('Mz', 'aveop1(beam.RMz)');

%ux=0;uy=0;
model.physics('beam').create('pdr01', 'DispRot0', 0);
model.physics('beam').feature('pdr01').selection.set([NodeNum]);
model.physics('beam').feature('pdr01').setIndex('Direction', true, 0);
model.physics('beam').feature('pdr01').setIndex('U0', 'ux', 0);
model.physics('beam').feature('pdr01').setIndex('Direction', true, 1);
model.physics('beam').feature('pdr01').setIndex('U0', 'uy', 1);
model.physics('beam').feature('pdr01').setIndex('DirectionRot', true, 0);


% fpout=fopen(datFile,'w+');
% fprintf(fpout,'%12s\t%12s\t%12s\t%12s\t%12s\t%12s\t','dis','Fx','Fy','Mz','F_Tang','F_Perp');
% fprintf(fpout,'\n');

Disp = 0.008;
theta = -pi/2;
Delta_F = 0;
Load = Disp*[0.5 0.6 0.65 0.7 0.75 0.8 0.85 1.0];
LoadN = length(Load);
F = zeros(LoadN,1);
for i=1:LoadN
    u_x=Load(i)*cos(theta);
    u_y=Load(i)*sin(theta);
    model.param.set('ux',u_x);
    model.param.set('uy',u_y);
    
    model.sol('sol1').runAll;
    
    Fx=mphglobal(model,'Fx');
    Fy=mphglobal(model,'Fy');
    %Mz=mphglobal(model,'Mz');
    %沿位移方向合力和垂直位移方向合力
    %F_Tang=Fx*cos(theta)+Fy*cos(pi/2-theta);
    %F_Perp=-Fx*sin(theta)+Fy*sin(pi/2-theta);
    Transf_Matrix=[cos(theta) sin(theta);
                   -sin(theta) cos(theta)];
    FF=Transf_Matrix*[Fx Fy]';
    F(i)=FF(1);
    %Delta_F=Delta_F+abs(Target_F(i)-F(i))/Target_F(i);
end
if F(1)>F(2)
    Delta_F=Delta_F+100*abs(Target_F(1)-F(1))/Target_F(1);
else
    Delta_F=Delta_F+abs(Target_F(1)-F(1))/Target_F(1);
end

for i=2:LoadN-1
    if F(i)<F(i-1)||F(i)>F(i+1)
        Delta_F=Delta_F+100*abs(Target_F(i)-F(i))/Target_F(i);
    else
        Delta_F=Delta_F+abs(Target_F(i)-F(i))/Target_F(i);
    end
end

if F(LoadN)<F(LoadN-1)
    Delta_F=Delta_F+100*abs(Target_F(LoadN)-F(LoadN))/Target_F(LoadN);
else
    Delta_F=Delta_F+abs(Target_F(LoadN)-F(LoadN))/Target_F(LoadN);
end
%Delta_F=Delta_F+2*(max(F)-min(F));
% fclose(fpout);
%Delta_F=Delta_F/0.05;
end
disp(Delta_F)
disp([num2str(toc),'  FINISHING CALCULATION'])
end