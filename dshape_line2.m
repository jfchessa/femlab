function dNdxi=dshape_line2(xi)

% function dNxi=dshape_line2(xi)
%
% Computes the gradient of the shape function dNxi for a 2 node line 
% element with respect to the parent coordinate system (xi)
%
%           1--------------------2
%        xi = -1              xi = 1
%
%    xi - the coordinate in the parent element space 
%
% function dNxi=dshape_line2() computes dNxi at the centroid
%
%
% Written by Jack Chessa, jfchessa@utep.edu

dNdxi=[-0.5;0.5];