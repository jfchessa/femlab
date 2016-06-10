function dNxi=dshape_tria6(coord)

% function dNxi=dshape_tria6(xi)
%
% Computes the gradient of the shape function,
% with respect to the parent coordinate system, DNxi
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
%    xi - the coordinate in the parent element space 
%
% function dNxi=dshape_tria6()
%
% Computes at the centriod
%
% Written by Jack Chessa, jfchessa@utep.edu

xi=coord(1);
eta=coord(2);
dNxi=[4*(xi+eta)-3   4*(xi+eta)-3;
             4*xi-1              0; 
                  0        4*eta-1;
     4*(1-eta-2*xi)          -4*xi;
              4*eta           4*xi;
             -4*eta  4*(1-xi-2*eta)];