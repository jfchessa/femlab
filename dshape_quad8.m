function dNdxi=dshape_quad8(coord)

% function dNdxi=dshape_quad8(xi)
%
% Computes the shape function N for a 4 node quadrilateral element. 
%
%    
%  (-1,1)    4---------7----------3 (1,1) 
%            |                    |
%            |                    |
%            |                    |
%            |                    |
%            8                    6
%            |                    |
%            |                    |
%            |                    |
%            |                    |
%  (-1,-1)   1---------5----------2%   (1,-1) 
%
%    xi - the coordinate in the parent element space to comopute N at
%
% function dNdxi=dshape_quad8() computes N at the centroid
%
%
% Written by Jack Chessa, jfchessa@utep.edu

if ( nargin==0 )
    s=0; t=0;
else
    s=coord(1); t=coord(2);
end
 
dNdxi = 0.25*[ -2*s*t - t^2 + 2*s + t ,  -s^2 - 2*s*t + s + 2*t ;
-2*s*t + t^2 + 2*s - t ,  -s^2 + 2*s*t - s + 2*t ;
2*s*t + t^2 + 2*s + t ,  s^2 + 2*s*t + s + 2*t ;
2*s*t - t^2 + 2*s - t ,  s^2 - 2*s*t - s + 2*t ;
4*s*t - 4*s ,  2*s^2 - 2 ;
-2*t^2 + 2 ,  -4*s*t - 4*t ;
-4*s*t - 4*s ,  -2*s^2 + 2 ;
2*t^2 - 2 ,  4*s*t - 4*t ];

