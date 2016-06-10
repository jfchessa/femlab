function N=shape_quad4(coord)

% function N=shape_quad4(xi)
%
% Computes the shape function N for a 4 node quadrilateral element. 
%
%    
%  (-1,1)    4--------------------3 (1,1) 
%            |                    |
%            |                    |
%            |                    |
%            |                    |
%            |                    |
%            |                    |
%            |                    |
%            |                    |
%            |                    |
%  (-1,-1)   1--------------------2%   (1,-1) 
%
%    xi - the coordinate in the parent element space to comopute N at
%
% function N=shape_quad4() computes N at the centroid
%
%
% Written by Jack Chessa, jfchessa@utep.edu

if ( nargin==0 )
    xi=0; eta=0;
else
    xi=coord(1); eta=coord(2);
end
      N=0.25*[ (1-xi)*(1-eta);
              (1+xi)*(1-eta);
              (1+xi)*(1+eta);
              (1-xi)*(1+eta)];