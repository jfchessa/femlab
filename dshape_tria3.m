function dNxi=dshape_tria3(coord)

% function dNxi=dshape_tria3(xi)
%
% Computes the gradient of the shape function,
% with respect to the parent coordinate system, DNxi
%
%    3  (0,1)
%    |\
%    | \
%    |  \
%    |   \
%    |    \
%    |     \
%    |      \
%    |       \
%    1--------2  (1,0)
%   (0,0) 
%
%    xi - the coordinate in the parent element space 
%
% function dNxi=dshape_tria3()
%
% Computes at the centriod
%
% Written by Jack Chessa, jfchessa@utep.edu

dNxi=[ -1 -1;1 0;0 1];