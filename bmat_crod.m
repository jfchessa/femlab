function B = bmat_crod( coord, c )

% B = bmat_crod( coord, c=1 )
%
% Generates B matrix of a 3D truss element with 
% torsional stiffness (similar to a NASTRAN CROD element)
%
% {      axial strain;
%  torsion eng shear strain}   = [B]*{d}
% 
% 	coord = coordinates at the element ends in column format.  The first column is the first 
%           node and second column holds the second node.
%   c = (optional, default=1) torsional stress recovery point.  If left unassigned you will 
%       need to multiply the second value in the strain to recover the correct torsinoal 
%       shear strain.

if nargin < 2
    c = 1;
end

invL = 1/norm(coord(:,1)-coord(:,2));
nx = (coord(1,2)-coord(1,1))*invL;
ny = (coord(2,2)-coord(2,1))*invL; 
nz = (coord(3,2)-coord(3,1))*invL;

n = [nx, ny, nz]*invL;

B = zeros(2,12);

B(1,1:3) = -n; 
B(1,7:9) =  n;
B(2,4:6) =  -c*n; 
B(2,10:12) = c*n;


       