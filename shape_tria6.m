function N=shape_tria3(xi)

% function N=shape_tria3(xi)
%
% Computes the shape function N for a 3 node triangular element. 
%
%    
%    3  (0,1)
%    |\
%    | \
%    |  \
%    |   \
%    6    5
%    |     \
%    |      \
%    |       \
%    1---4----2  (1,0)
%   (0,0) 
%
%    xi - the coordinate in the parent element space to comopute N at
%
% function N=shape_tria6() computes N at the centroid
%
%
% Written by Jack Chessa, jfchessa@utep.edu

if ( nargin==0 )
    xi=[0.3333333333333333333,0.33333333333333333333333];
end

N=[1-3*(xi(1)+xi(2))+4*xi(1)*xi(2)+2*(xi(1)*xi(1)+xi(2)*xi(2));
    xi(1)*(2*xi(1)-1);
    xi(2)*(2*xi(2)-1);
    4*xi(1)*(1-xi(1)-xi(2));
    4*xi(1)*xi(2);
    4*xi(2)*(1-xi(1)-xi(2)) ];