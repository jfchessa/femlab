function order = element_porder( etype )

% function ORDER = ELEMENT_PORDER( ETYPE )
%
% Returns the interpolation polynomial order of an element's shape
% function.
%
% THis is part of FEMLAB
% Written by Jack Chessa, jfchessa@utep.edu

switch lower(etype)
    
    case {'point1','none','unknown'}
        order=0;
    
    case {'line2','tria3','quad4','tetra4','hexa8','prism6','pyramid5'}
        order=1;
        
    case {'line3','tria6','quad8','quad9','tetra10','hexa20','hexa27',...
            'prism18','prism15','pyramid14','pyramid13'}
        order=2;
           
    otherwise
        order=0;
end

end

