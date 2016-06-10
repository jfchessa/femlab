function [ke,A] = kmat_tria3(coord, E, nu, thk)

%  [ke,A] = kmat_tria3(coord, E, nu, thk)
% Generates equations for a plane stress tria3 element, ke 
% as well as the area of the elemetn A.
% E = Youngs modulus
% nu = Poissons ratio
% thk = element thickness
% coord = coordinates at the element ends

c1=E/(1-nu^2);  % plane stress material stiffness matrix
c2=nu*c1;
c3=0.5*(1-nu)*c1;
C=[c1 c2 0;c2 c1 0; 0 0 c3];
[B,A]=bmat_tria3(coord);
ke=B'*C*B*A*thk;
