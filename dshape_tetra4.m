function dNxi=dshape_tetra4(xi)

% function N=dshape_tetra4(xi)
%
% Computes the derivative of the shape function dNxi with respect to the
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
% function N=dshape_tetra4() computes N at the centroid
%
%
% Written by Jack Chessa, jfchessa@utep.edu


dNxi=[ -1 -1 -1;1 0 0;0 1 0; 0 0 1 ];