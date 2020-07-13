function [x,f1]=FR_ConjGrad(F,gradF,x0,TolGrad,MaxIter,h,n,ep)
%
%  Input
%  F........ function handle for objective function F(x) with input argument, x
%  gradF... function handle for gradient of objective function with input argument, x
%  x0....... design variable vector initial guess for starting point
%  TolGrad.. Tolerance for norm of objective gradient to be zero
%  MaxIter.. Maximum number of iterations
%  h.........f'(x)=(f(x+h)-f(x-h))/(2*h)
%  n.........nth element of x
%  ep........ep is used for adjusting the step lenght alpha.
%
%  Output
%  x........ final design variable vector found to minimize objective function, F(x)
%  f1....... final objective function minimum value

if nargin<5 || isempty(TolGrad), TolGrad=1e-2; end
if nargin<4 || isempty(MaxIter), MaxIter=20; end

%% Change following steepest descent algorithm to Conjugate Gradient Method

% Initialize loop parameters
iter = 0;
f0   = F(x0);
c    = gradF(x0,h,n);
Converged = norm(c) < TolGrad;
disp('   iter  alpha     f(x)         norm(c)    x[n]')
fprintf(' %6.0f %5.3f %8.4f %8.4f %5.4f\n',[iter, 0, f0, norm(c),x0(n)])
alpha = ep*norm(c);
f_old=f0;

%% Search direction and line search iterations
   d=-c;
   x_old=x0;
   f = @(alpha) F(x0+alpha*d);
   alphaUpper = bracket( f, 0, 0.1*alpha );
   [alpha,f1] = fminbnd( f, 0, alphaUpper );
   x = x0 + alpha*d;
   c_before=c;
   c = gradF(x0,h,n);
   Converged = (norm(c) < TolGrad);
   x0 = x;
   d_before=d;
   iter=iter+1;
      if(f_old<f1)
       f1=f_old;
       x=x_old;
       return
   end
   if(f1==-1)
       x=x_old;
       return
   end
   f_old=f1;
   x_old=x;
   fprintf(' %6.0f %5.3f %8.4f %8.4f %5.4f\n',[iter, 0, f1, norm(c),x0(n)])
while iter<MaxIter && ~Converged
   iter = iter + 1;
   beta=(c'*(c-c_before))/(c_before'*c_before);
   %beta=(c'*c)/(c_before'*c_before);
   d=-c+beta*d_before;
   f = @(alpha) F(x0+alpha*d);
   alphaUpper = bracket( f, 0, 0.1*alpha );
   [alpha,f1] = fminbnd( f, 0, alphaUpper );
   
   if(f1==-1)
       f1=f_old;
       x=x_old;
       break;
   end  %如果超出范围则退出
   if(f1>f_old)
       x=x_old;
       break;
   end %大于上次则退出
   x = x0 + alpha*d;
   c_before=c;
   c = gradF(x0,h,n);
   Converged = (abs(f1-f_old)<0.001);
   x0 = x;
   d_before=d;   
   f_old=f1;
   x_old=x;

   fprintf(' %6.0f %5.3f %8.4f %8.4f %5.4f\n',[iter, alpha, f1, norm(c),x(n)])
end
end

%%Bracket interval for 1-D line search
function [alphaUpper,alphaLower,nfunc] = bracket( f, alpha0, delta, MaxIter )
	% usage: [alphaUpper,alphaLower,nfunc] = bracket( f, alpha0, delta, MaxIter )
	%  Golden section search to bracket unimodal, univariate minimum
	%--input
	%  f = function handle to univariate function of alpha
	%  alpha0... Starting point lower bound on bracket
	%  delta.... Guess for upper bound on bracket on unimodal min
	%  MaxIter.. Maximum number of iterations
	%--output
	%  alphaUpper... Upper bound on alpha to bracket min of f(alpha)
	%  alphaLower... Lower bound on alpha to bracket min of f(alpha)
	%  nfunc........ Number of function evaluations

      if nargin<3, delta = 1.e-3; end
      if nargin<4, MaxIter=10; end
      %--Local variables
      %  expand... Expansion factor for extending Upper Bound
      %  phi...... Golden section ratio = 1.681...
      phi=(1+sqrt(5))/2;
      expand = phi; % Set=1 to use equal interval search
      % Initialize variables
      alphaLower=alpha0;
      alphaLast=alpha0;
      alphaNext=delta;
      fLast=f(alphaLast);
      fNext=f(alphaNext);
      if(fNext==-1)
          alphaUpper=alpha0;
      end
      iter=1;
      while fNext<fLast && iter<=MaxIter
         iter = iter + 1;
         delta = expand*delta;
         alphaLower = alphaLast;
         alphaLast = alphaNext;
         alphaNext = alphaNext + delta;
         fLast = fNext;
         fNext = f(alphaNext);
         if (fNext==-1)
             alphaNext=alphaLast;
             break;
         end
      end
      alphaUpper=alphaNext;
      nfunc=iter+1;
   end