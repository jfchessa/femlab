function N=shape_hexa8(coord)

% function N=shape_hexa8(xi)
%
% Computes the shape function N for a 8 node hexahedral element. 
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
% function N=shape_hexa8() computes N at the centroid
%
%
% Written by Jack Chessa, jfchessa@utep.edu

if ( nargin==0 )   
    I1=1/2;
    I2=1/2;

else
    I1=1/2-coord/2;
    I2=1/2+coord/2;
end
     
      N=[   I1(1)*I1(2)*I1(3);
            I2(1)*I1(2)*I1(3);
            I2(1)*I2(2)*I1(3);
            I1(1)*I2(2)*I1(3);
            I1(1)*I1(2)*I2(3);
            I2(1)*I1(2)*I2(3);
            I2(1)*I2(2)*I2(3);
            I1(1)*I2(2)*I2(3)   ];