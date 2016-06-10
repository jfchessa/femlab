function N=shape_tetra4(xi)

% function N=shape_tetra4(xi)
%
% Computes the  shape function N with respect to the
% parent coordinate system (xi,eta,zeta) for a 3 node triangular element. 
%
%    
%                4   (0,0,1)
%              / | \
%             /  |  \
%            /   |   \ 
%           /    |    \ 
%          /     |     \
% (0,0,0) 1 -----|------3 (0,1,0)
%            -   2  -
%               (1,0,0)
%              
%
%    xi - the coordinate in the parent element space to comopute N at
%
% function N=shape_tetra4() computes N at the centroid
%
%
% Written by Jack Chessa, jfchessa@utep.edu


if ( nargin==0 )
    xi=[ 0.166666666666666666667, 0.166666666666666666667, 0.166666666666666666667];
end

N=[ 1-xi(1)-xi(2)-xi(3); xi(1); xi(2); xi(3) ];