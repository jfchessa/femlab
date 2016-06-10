function d=absr(x,eps)

% function  d=absr(x,eps)
%
% returns a regularized absolute value function of the variable
% x.  The discontinutiy at x=0 is smoothed over over x=[-eps,eps]

d=abs(x);

is=find(d<=eps);

a0=1/16/eps*(15-14*eps+2*eps^2);
a2=-3/8/eps^2*(5-10*eps+2*eps^2);
a4=5/16/eps^5*(3-6*eps+2*eps^2);

d(is)=a0+a2*x(is).^2+a4*x(is).^4;

% plot(x,a0+a2*x.^2+a4*x.^4)
% hold on
% plot([-2 0 2],[2 0 2],'g')