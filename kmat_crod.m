function ke = kmat_crod( coord, AE, JG )

% k = kmat_crod( coord, AE, JG )
%
% Generates stiffness matrix of a 3D truss element with 
% torsional stiffness (similar to a NASTRAN CROD element)
% 	AE = modulus of elasticity
% 	JG = Area of cross-section
% 	coord = coordinates at the element ends

x1=coord(1,1); y1=coord(1,2); z1=coord(1,3);
x2=coord(2,1); y2=coord(2,2); z2=coord(2,3);

L=sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);
nx=(x2-x1)/L; ny=(y2-y1)/L; nz=(z2-z1)/L;

n = [nx, ny nz];
N = n'*n;
Z = zeros(3,3);

ka = (AE/L)*N;
kt = (JG/L)*N;

ke = [  ka    Z  -ka    Z;
         Z   kt    Z  -kt;
	   -ka    Z   ka    Z;
         Z  -kt    Z   kt ];
       