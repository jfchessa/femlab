function d=heavisider(eta,m,k)

% function  d=heavisider(eta,m,k)
%
% returns a regularized Heaviside step function of the variable
% x.  The true Heaviside function is approximated by a function
% which is has a smooth transition over the open set (-eps,eps).
%
% k is the continuity of the regularization and m is the number
% of vanishing moments.  The default is m=2 k=1. Other valid
% options are {m,k}={0,0},{2,0},{2,1},{2,2},{4,0},{4,1},{4,2}

d=heaviside(eta);

d(find(eta>1))=1;

is=find(abs(eta)<1);

switch 10*m+k
  
case 00
  d=heaviside(eta);
  
case 20
  d(is)=1/2+1/8*(9*eta(is)-5*eta(is).^3);
  
case 22
  d(is)=1/2+1/64*(105*eta(is)-175*eta(is).^3+147*eta(is).^5-45*eta(is).^7);
  
case 40
  d(is)=1/2+1/128*(225*eta(is)-350*eta(is).^3+189*eta(is).^5);
  
case 41
  d(is)=1/2+1/256*(525*eta(is)-1225*eta(is).^3+...
                                1323*eta(is).^5-495*eta(is).^7);
  
case 42
  d(is)=1/2+1/2048*(4725*eta(is)-14700*eta(is).^3+23814*eta(is).^5 ...
                             -17820*eta(is).^7+5005*eta(is).^9);
  
otherwise  % 2 1
  d(is)=1/2+1/32*(45*eta(is)-50*eta(is).^3+21*eta(is).^5);
  
end
  