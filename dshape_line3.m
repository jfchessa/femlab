function dNdxi=dshape_line3(xi)

% function dNxi=dshape_line3(xi)
%
% Computes the gradient of the shape function dNxi for a 3 node line 
% element with respect to the parent coordinate system (xi)
%
%           1---------3----------2
%        xi = -1              xi = 1
%
%    xi - the coordinate in the parent element space 
%
% function dNxi=dshape_line3() computes dNxi at the centroid
%
%
% Written by Jack Chessa, jfchessa@utep.edu

if ( nargin==0 )
    xi=0.0;
end

dNdxi=[ xi-.5; -2*xi; xi+.5 ];