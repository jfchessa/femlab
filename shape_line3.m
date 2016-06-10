function N=shape_line3(xi)

% function N=shape_line3(xi)
%
% Computes the shape function N for a 3 node line element. 
%
%           1---------3----------2
%        xi = -1              xi = 1
%
%    xi - the coordinate in the parent element space to comopute N at
%
% function N=shape_line3() computes N at the centroid
%
%
% Written by Jack Chessa, jfchessa@utep.edu

if ( nargin==0 )
    xi=[0.0];
end

N=[ -0.5*(1-xi)*xi; 1-xi^2; 0.5*(1+xi)*xi ];