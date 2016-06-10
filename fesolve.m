function [d,freac]=fesolve(K,f,ifix,ival)

% function [d,freac]=fesolve(K,f,ifix,ival)
%
% Solve a linear system Kd=f with the constraints that d(ifix)=ival.
% freac is the constraint forces.
%
% function [d,freac]=fesolve(K,f,ifix)
%
% Solve a linear system Kd=f with the constraints that d(ifix)=0.
% freac is the constraint forces.
%
% Written by Jack Chessa jfchessa@utep.edu
%

nc=length(ifix);  % number of dof constraints
nn=length(f);     % number of dofs

if ( nc==0 )
    d=K\f;
    freac=[];
    return
end

isolve=setdiff(1:nn,ifix);

d=zeros(nn,1);
freac=zeros(nc,1);

if ( nargin==4 ) % then we have inhomogenous bc
  if ( size(ival,1) == 1 ) % ival is a row vector
    f = f - K(:,ifix)*ival';
    d(ifix)=ival';
  else % ival is a column vector
    f = f - K(:,ifix)*ival;
    d(ifix)=ival;
  end
  
end

d(isolve)=K(isolve,isolve)\f(isolve);  % solve

f=K*d;
freac=f(ifix);