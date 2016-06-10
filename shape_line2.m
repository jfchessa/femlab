function N=shape_line2(xi)

% function N=shape_line2(xi)
%
% Computes the shape function N for a 2 node line element. 
%
%           1--------------------2
%        xi = -1              xi = 1
%
%    xi - the coordinate in the parent element space to comopute N at
%
% function N=shape_line2() computes N at the centroid
%
%
% Written by Jack Chessa, jfchessa@utep.edu

if ( nargin==0 )
    xi=[0.0];
end

N=0.5*[ 1-xi; 1+xi ];