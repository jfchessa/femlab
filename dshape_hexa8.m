function dNdxi=dshape_hexa4(coord)

% function dNdxi=dshape_hexa4(xi)
%
% Computes the gradient of the shape function dNdxi for a 8 node hexahedral 
% element. 
%
%    
% 
%                        8 (-1,1,1)
%                     /    \    
%                  /          \
%               /                \
%  (-1,-1,1) 5                     \
%            |\                     7 (1,1,1)
%            |   \    (-1,1,-1)   / |
%            |     \     4    /     |
%            |        \    /        |
%            |           6 (1,-1,1) | 
% (-1,-1,-1) 1           |          |
%             \          |          3(1,1,-1)
%                \       |        /
%                  \     |     /
%                     \  |  /
%                        2(1,-1,-1)
%
%
%    xi - the coordinate in the parent element space to comopute N at
%
% function dNdxi=dshape_hexa4() computes N at the centroid
%
%
% Written by Jack Chessa, jfchessa@utep.edu

if ( nargin==0 )
    xi=0; eta=0; zeta=0;   

else
    xi=coord(1); eta=coord(2); zeta=coord(3);
end
     
     dNdxi=[   -1+eta+zeta-eta*zeta   -1+xi+zeta-xi*zeta  -1+xi+eta-xi*eta;
                 1-eta-zeta+eta*zeta   -1-xi+zeta+xi*zeta  -1-xi+eta+xi*eta;
                 1+eta-zeta-eta*zeta    1+xi-zeta-xi*zeta  -1-xi-eta-xi*eta;
                -1-eta+zeta+eta*zeta    1-xi-zeta+xi*zeta  -1+xi-eta+xi*eta;      
                -1+eta-zeta+eta*zeta   -1+xi-zeta+xi*zeta   1-xi-eta+xi*eta;
                 1-eta+zeta-eta*zeta   -1-xi-zeta-xi*zeta   1+xi-eta-xi*eta;
                 1+eta+zeta+eta*zeta    1+xi+zeta+xi*zeta   1+xi+eta+xi*eta;
                -1-eta-zeta-eta*zeta    1-xi+zeta-xi*zeta   1-xi+eta-xi*eta  ]/8;