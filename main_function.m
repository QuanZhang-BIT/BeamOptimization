clc;clear;
tic;

A=[];
b=[];
Aeq=[];
beq=[];
IntCon=9;
lb=[0 0 0 0 0 0 0 0  10  0.1 0.1 0.1 0.1 0.1 0.1];
ub=[0.03 0.015 0.03 0.015 0.03 0.015 0.03 0.015  10  5.0 5.0 5.0 5.0 5.0 5.0];
rng default
fun = @ResidualF;
nonlcon = @constraint;
options = optimoptions('ga','PopulationSize',200,'MaxGenerations',100,'PlotFcn', @gaplotbestf);
[x,fval,exitflag,output,population,scores] = ga(fun,15,A,b,Aeq,beq,lb,ub,nonlcon,IntCon,options)

disp([num2str(toc),'  FINISHING CALCULATION'])

