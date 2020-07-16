%x(1) -0.2339  33.9404
%x(2) -0.1513  0.2439
%x(3) 0.0702   16.2695
%x(4) 5e8      inf
%x(5) 1e-5
%x(6) inf
%x(7) inf

count=0;
a=1e8;
b=4e9;
x = [14.78;1.2e-1;1.5e-1;4e9;1e-5;4e-3;5e-12;9.5];
element=2;
temp=-0.15:0.01:0.24;
y=temp;
count=1;
for i=temp
    x(element)=i;
    y(count)=equationfft(x);
    count=count+1;
end

plot(temp,y);
% 
% while(abs(a-b)>1e6)
%     x(element)=a;
%     re_a=equationfft(x);
%     x(element)=b;
%     re_b=equationfft(x);
%     if (re_a==-1 && re_b~=-1)||(re_a~=-1 && re_b==-1)
%         mid=(a+b)/2;
%         x(element)=mid;
%         re_m=equationfft(x);
%         if(re_m~=-1)
%             a=mid;
%         else 
%             b=mid;
%         end
%     else
%         disp('error');
%     end
% end

    
        
        