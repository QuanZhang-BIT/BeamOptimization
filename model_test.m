clc;clear;
tic;
model=mphload('model_test.mph');
x=[0.002230521	 0.005082751	 0.012339819	 0.003589266	 0.018917285	 0.008258576	 0.029156795	 0.004054656];
in=length(x);
in=in/2;
%control points
P=zeros(2,in+2);
P(1:2,1)=(0);%起点
P_end = [0.03;0.015];
P(1:2,in+2)=P_end;%终点

for i=1:in
    P(1,i+1)=x(2*i-1);%中间控制点x坐标
    P(2,i+1)=x(2*i);%中间控制点y坐标
end
n = 3; %curve's degree n = p + 1 where p is the degree of the curve
%knot vector
t = [0 0 0 0.3 0.5 0.7 1 1 1];
%weight vector
w = [1.579418929	 2.379914318	 1.701708739	 3.086466483	 3.952908474	 3.371585436];
%w = [1 1 1 1 1 1];
%call of the function
C = nurbsfun(n, t, w, P);

nOrder=size(C,2);
for ii=1:nOrder
    strX_value=C(1,ii);strY_value=C(2,ii);
    model.geom('geom1').feature('pol1').setIndex('table', strX_value, ii-1, 0);
    model.geom('geom1').feature('pol1').setIndex('table', strY_value, ii-1, 1);
end
model.geom('geom1').run('pol1');
model.geom('geom1').run;

model.param.set('b',0.001);
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

datFile='0-0_3.0-1.5cm(disp0.8(0.5-1.0))temp.txt';
fpout=fopen(datFile,'w+');
fprintf(fpout,'%12s\t%12s\t%12s\t%12s\t%12s\t%12s\t','dis','Fx','Fy','Mz','F_Tang','F_Perp');
fprintf(fpout,'\n');

Disp = 0.015;
theta = -pi/2;
Delta_F = 0;
Load = Disp*[0.005:0.005:0.1 0.11:0.01:0.13 0.14:0.02:1.0];
%Load = Disp*[0.6 0.65 0.7 0.75 0.8 0.85 1.0];
LoadN = length(Load);
for i=1:LoadN
    u_x=Load(i)*cos(theta);
    u_y=Load(i)*sin(theta);
    model.param.set('ux',u_x);
    model.param.set('uy',u_y);
    
    model.sol('sol1').runAll;
    
    Fx=mphglobal(model,'Fx');
    Fy=mphglobal(model,'Fy');
    Mz=mphglobal(model,'Mz');
    %沿位移方向合力和垂直位移方向合力
    %F_Tang=Fx*cos(theta)+Fy*cos(pi/2-theta);
    %F_Perp=-Fx*sin(theta)+Fy*sin(pi/2-theta);
    Transf_Matrix=[cos(theta) sin(theta);
                   -sin(theta) cos(theta)];
    F=Transf_Matrix*[Fx Fy]';
    fprintf(fpout,'%12.6f\t%12.6f\t%12.6f\t%12.6f\t%12.6f\t%12.6f\t',Load(i),Fx,Fy,Mz,F(1),F(2));
    fprintf(fpout,'\n');
end
fclose(fpout);
mphsave(model,'Temp')
disp([num2str(toc),'  FINISHING CALCULATION'])

hold on
plot(C(1, :), C(2, :), 'Color',[0.9883 0.3438 0.1445],'LineWidth',6);
plot(P(1, 2:end-1), P(2, 2:end-1), 'kx','LineWidth',3,'MarkerSize',12);
%plot(P(1, 2:end), P(2, 2:end), 'kx','LineWidth',3,'MarkerSize',12);
plot(P(1, :), P(2, :), 'k','LineWidth',1);
%plot(C(1, :), C(2, :), 'k','LineWidth',2);

set(gca,'FontSize',16);
axis equal;
axis off;

