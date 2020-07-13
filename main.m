format long g  %取消科学计数法
% CHO = 4e-1;
% CHCO3 = 1.2e-1;
% CCO3 = 1.5e-1;
% K1 = 4e9;
% ka = 1e-5;
% kC = 4e-3;
% kp = 5e-12;
% p = 9.5;
%目前已经得出的两组数据
%result<0.1
%x=[-0.1;-0.01;0.15;4000000000.00000;1.0000000000000e-05;0.00400000000000000;5.00000000000000e-12;9.50000000000000]
%x=[18.2183246425871;-0.0502910181963702;0.15;4000000000;1e-05;0.004;5e-12;9.5];
x = [4e-1;0.12;1.5e-1;4e9;1e-5;4e-3;5e-12;9.5];
%x_con=[-36.5;-0.0327;-5.86e2;4.1e9;1e-5;-1.18e3;-9.3e12;9.5];
%x_th=[0.41;0.13;0.16;4.1e9;1e-5;4.1e-3;5.1e-12;9.5];
%
%x=[-0.327715618592705;1.2e-1;1.5e-1;4e9;1e-5;4e-3;5e-12;9.5];---1
%x=[-0.327715618592705;0.019;1.5e-1;4e9;1e-5;4e-3;5e-12;9.5];---2
%x=[-0.327715618592705;0.019;0.147;4e9;1e-5;4e-3;5e-12;9.5];---3
%x=[-0.327715618592705;0.019;0.147;1041181460.061;1e-5;4e-3;5e-12;9.5];---4
%-----------------------------------------------------------------------5
%
%

result=100;
while result>4
    x_old=x;
for n=1:7
    switch(n)
    case 1 
     fprintf('%dth element\n',n);
    [x_new,result]=FR_ConjGrad(@(x) equationfft(x),@(x,h,n) diff_by_h(x,h,n),x,1e-3,20,0.01,n,1e-11);  
    fprintf('%f --> %f\n',[x(n),x_new(n)]);
    x=x_new;
    case 2 
     fprintf('%dth element\n',n);
    [x_new,result]=FR_ConjGrad(@(x) equationfft(x),@(x,h,n) diff_by_h(x,h,n),x,1e-3,20,0.01,n,1e-18);
    fprintf('%f --> %f\n',[x(n),x_new(n)])
    x=x_new;
    end
    disp(x);
end
delta=abs(x-x_old);

if delta(1)<0.05
    x(1)=x(1)+(rand-0.5)*0.1;
    if x(1)>33.9404 || x(1)<-0.2339
        x(1)=0.4;
    end
end
if delta(2)<0.05
    x(2)=x(2)+(rand-0.5)*0.2;
    if x(2)<-0.1513 || x(2)>0.2439
        x(2)=0.12;
    end
end
end

for n=3:7
    switch n
    case 3
    fprintf('%dth element\n',n);
    [x_new,result]=FR_ConjGrad(@(x) equationfft(x),@(x,h,n) diff_by_h(x,h,n),x,1e-3,5,0.01,n,1e-21);
    fprintf('%f --> %f\n',[x(n),x_new(n)])
    x=x_new;
    case 4
    fprintf('%dth element\n',n);
    [x_new,result]=FR_ConjGrad(@(x) equationfft(x),@(x,h,n) diff_by_h(x,h,n),x,1e-3,20,0.0001,n,1e16);
    fprintf('%f --> %f\n',[x(n),x_new(n)])
    x=x_new;
    case 5
    continue;
    case 6
    fprintf('%dth element\n',n);
    [x_new,result]=FR_ConjGrad(@(x) equationfft(x),@(x,h,n) diff_by_h(x,h,n),x,1e-3,10,0.01,n,1e-16);
    fprintf('%f --> %f\n',[x(n),x_new(n)])
    x=x_new;
    case 7
    fprintf('%dth element\n',n);
    [x_new,result]=FR_ConjGrad(@(x) equationfft(x),@(x,h,n) diff_by_h(x,h,n),x,1e-3,20,0.01,n,1e-32);
    fprintf('%f --> %f\n',[x(n),x_new(n)])
    x=x_new;
    end
end