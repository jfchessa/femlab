function [s,nne]=etypestr(etype)

% function [ETYPE,NNE]=ETYPESTR(E)
%
% Converts an element type id as from GMSH to an actual element type string
% NNE - returns teh number of nodes per element

switch etype
    case 1
        s='Line2';
        nne=2;
    case 2
        s='Tria3';
        nne=3;
    case 3
        s='Quad4';
        nne=4;
    case 4
        s='Tetra4';
        nne=4;
    case 5
        s='Hexa8';
        nne=8;
    case 6
        s='Prism6';
        nne=6;
    case 7
        s='Pyramid5';
        nne=5;
    case 8
        s='Line3';
        nne=3;
    case 9
        s='Tria6';
        nne=6;
    case 10
        s='Quad9';
        nne=9;
    case 11
        s='Tetra10';
        nne=10;
    case 12
        s='Hexa27';
        nne=27;
    case 13
        s='Prism18';
        nne=18;
    case 14
        s='Pyramid14';
        nne=14;
    case 15
        s='Point1';
        nne=1;
    case 16
        s='Quad8';
        nne=8;
    case 17
        s='Hexa20';
        nne=20;
    case 18
        s='Prism15';
        nne=15;
    case 19
        s='Pryamid13';
        nne=13;
    case 20
        s='TriInc9';
        nne=9;
    case 21
        s='Tria10';
        nne=10;
    case 22
        s='TriaInc12';
        nne=13;
    case 23
        s='Tria15';
        nne=15;
    case 24
        s='TriaInc15';
        nne=15;
    case 25
        s='Tria21';
        nne=21;
    case 26
        s='Line4';
        nne=4;
    case 27
        s='Line5';
        nne=5;
    case 28
        s='Line6';
        nne=6;
    case 29
        s='Tetra20';
        nne=20;
    case 30
        s='Tetra35';
        nne=35;
    case 31
        s='Tetra56';
        nne=56;
    otherwise
        s='Unknown';
        nne=0;
end
end
