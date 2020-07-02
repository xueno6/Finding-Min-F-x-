% CHO = 4e-1;
% CHCO3 = 1.2e-1;
% CCO3 = 1.5e-1;
% K1 = 4e9;
% ka = 1e-5;
% kC = 4e-3;
% kp = 5e-12;
% p = 9.5;
x = [4e-1;1.2e-1;1.5e-1;4e9;1e-5;4e-3;5e-12;9.5];
%初始值
%x(4)和x(8)无需改变

result=equationfft(x); 
%这是目标函数，求出 min equationfft(x)
%注意，若x超出equationfft(x)范围，则输出-1
