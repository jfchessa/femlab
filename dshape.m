function dNdxi=dshape(etype,coord)

% function DNDXI=DSHAPE(ETYPE,COORD)
%
% Computes the shape function N for various element types. 
%
%
%   ETYPE the element type. valid element types are :
%       Line2, Line3, Tria3, Tria6, Quad4, Quad8, Tetra4, Tetra10, 
%       Hexa8, Hexa27
%
%   COORD - the coordinate in the parent element space to comopute N at
%
% function DNDXI=DSHAPE(ETYPE) computes N at the centroid
%
%
% This function is part of FEMLAB
% Written by Jack Chessa, jfchessa@utep.edu


if ( nargin>1 )
    dNdxi=feval(['dshape_',lower(etype)],coord);
else
    dNdxi=feval(['dshape_',lower(etype)]);
end

end

