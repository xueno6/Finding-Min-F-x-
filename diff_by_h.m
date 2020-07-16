function grad=diff_by_h(x,h,n)
%f is the orginal function
%x is the point we need its diff
%h is the delta x
%n means the nth element of vector x
% we assume that x is a vector of 8*1
x_h=x;
x_h(n)=x_h(n)+h;
h_x=x;
h_x(n)=h_x(n)-h;
grad=zeros(8,1);
grad(n)=(equationfft(x_h)-equationfft(h_x))/(2*h);
end
