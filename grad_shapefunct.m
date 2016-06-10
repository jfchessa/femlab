function [dNx,jac]=grad_shapefunct(coord,dNdxi,sdim,edim)

%
% function [dNx,jac]=grad_shapefunct(coord,dNxi,sdim,edim)
%
% Computes the gradient of a shape function with respect to the spacial
% coordinate system, dNx and the Jacobian of the element, jac given
%
%    dNdxi -  gradient of the shape function w.r.t. parent cs
%    coord - nodal coordinate matrix (row oriented)
%    sdim - spacial dimension (default is # cols in coord)
%    edim - dimension of parent cs (default is #cols dNdxi)
%

if ( nargin<3 )
    sdim=size(coord,2);
end

if (nargin<4)
    edim=size(dNdxi,2);
end

[jac,jmat] = element_jacobian(dNdxi,coord,sdim,edim);

if ( sdim==edim )
    dNx = dNdxi*inv(jmat);
    
elseif ( edim==1 )
    invj=1./jmat';
    invj( find(invj==Inf) ) = 0;
    dNx = dNdxi*invj;
    
else % edim=2 and sdim=3
    x_xi = jmat(:,1);
    x_eta = jmat(:,2);
    jmat(:,3) = cross( x_xi, x_eta );
    jmat(:,3) = jmat(:,3)/norm(jmat(:,3));
    invj=inv(jmat);
    dNx = dNdxi*invj(1:edim,:);
    
end