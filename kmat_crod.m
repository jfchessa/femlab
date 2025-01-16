function ke = kmat_crod( coord, AE, JG )

% k = kmat_crod( coord, AE, JG )  
% Generates stiffness matrix of a 3D truss element with 
% torsional stiffness (similar to a NASTRAN CROD element)
% 	AE = modulus of elasticity
% 	JG = Area of cross-section
% 	coord = coordinates at the element ends in column format.  The first column is the first 
%           node and second column holds the second node.

invL = 1/norm(coord(:,1)-coord(:,2));
nx = (coord(1,2)-coord(1,1))*invL;
ny = (coord(2,2)-coord(2,1))*invL; 
nz = (coord(3,2)-coord(3,1))*invL;

n = [nx, ny nz];
N = n'*n;
Z = zeros(3,3);

ka = (AE*invL)*N;
kt = (JG*invL)*N;

ke = [  ka    Z  -ka    Z;
         Z   kt    Z  -kt;
	   -ka    Z   ka    Z;
         Z  -kt    Z   kt ];
       