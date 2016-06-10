function [d,freac]=fesolve2(K,f,ifix,ival)

% function [d,freac]=fesolve2(K,f,ifix,ival)
%
% Solve a linear system Kd=f with the constraints that d(ifix)=ival.
% freac is the constraint forces. Uses the equation substitution method
% presented in class.  (ifix and ival must be row vectors)
%
% function [d,freac]=fesolve2(K,f,ifix)
%
% Solve a linear system Kd=f with the constraints that d(ifix)=0.
% freac is the constraint forces.
%
% Written by Jack Chessa jfchessa@utep.edu
%

nc=length(ifix);  % number of dof constraints
nn=length(f);     % number of dofs

Kreac=K(ifix,:);  % store this part of the stiffness matrix to compute the
                  % reaction force.  We will directly modify K
                  
if ( nargin==4 ) % then inhomogenious bcs
  f=f-K(:,ifix)*ival';  % move columns of KI*dI to the rhs
end

K(:,ifix)=0.0;  % zero out columns
K(ifix,:)=0.0;  % zero out rows
K(ifix,ifix)=eye(nc); % put ones on the diagonal

f(ifix)=ival; % set rhs for dof constraint equations

d=K\f;  % solve for nodal dofs

freac=Kreac*d;  % compute reaction forces


