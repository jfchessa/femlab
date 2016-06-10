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
dNx = dNdxi*inv(jmat);