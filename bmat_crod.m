function B = bmat_crod( coord )

% B = bmat_crod( coord )
%
% Generates B matrix of a 3D truss element with 
% torsional stiffness (similar to a NASTRAN CROD element)
% 
% 	coord = coordinates at the element ends

x1=coord(1,1); y1=coord(1,2); z1=coord(1,3);
x2=coord(2,1); y2=coord(2,2); z2=coord(2,3);

L=sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);
nx=(x2-x1)/L; ny=(y2-y1)/L; nz=(z2-z1)/L;

n = [nx, ny nz]/L;

B = zeros(2,12);

B(1,1:3) = -n; 
B(1,7:9) =  n;
B(2,4:6) = -n; 
B(2,10:12) =  n;


       