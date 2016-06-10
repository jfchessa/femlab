function [ qpts, qwts ] = element_quadrature( etype, order )

% function [ QPTS, QWTS ] = ELEMENT_QUADRATURE( ETYPE, ORDER )
%
%     Returns the appropriate quadrature points and weigths for a
%     particular element type and the polynomial order of the integrand.
%
%   ETYPE - The element type string   
%   ORDER - The polynomial order of the integral
%
% This is part of FEMLAB
%
% written by Jack Chessa, jfchessa@utep.edu

%       1  Line2    2-node line.
%       2  Tria3    3-node triangle.
%       3  Quad4    4-node quadrangle.
%       4  Tetra4   4-node tetrahedron.
%       5  Hexa8    8-node hexahedron.
%       6  Prism6   6-node prism.
%       7  Pyramid5 5-node pyramid.
%       8  Line3    3-node second order line  
%       9  Tria6    6-node second order triangle  
%      10  Quad9    9-node second order quadrangle  
%      11  Tetra10  10-node second order tetrahedron  
%      12  Hexa27   27-node second order hexahedron  
%      13  Prism18  18-node second order prism  
%      14  Pyramid14 14-node second order pyramid  
%      15  Point1   1-node point.
%      16  Quad8    8-node second order quadrangle  
%      17  Hexa20   20-node second order hexahedron  
%      18  Prism15  15-node second order prism  
%      19  Pyramid13 13-node second order pyramid 

switch lower(etype)
   
    case {'line2','line3' }
        [ qpts, qwts ] = quadrature_gaussian( order, 1 ); 
        
    case {'quad4','quad8','quad9'}
        [ qpts, qwts ] = quadrature_gaussian( order, 2 );
        
    case {'hexa8','hexa20','hexa27'}
        [ qpts, qwts ] = quadrature_gaussian( order, 3 );
        
    case {'tria3','tria6'}
        [ qpts, qwts ] = quadrature_simplex( order, 2 );
        
    case {'tetra4','tetra10'}
        [ qpts, qwts ] = quadrature_simplex( order, 3 );    
        
    otherwise
        disp('Unsupported element type in ELEMENT_QUADRATURE')

end


end

