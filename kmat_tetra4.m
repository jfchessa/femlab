function ke = kmat_tetra4(coord, E, nu)

%  ke = kmat_tetra4(coord, E, nu)
% Generates equations for a tetra4 element, ke 
% E = Youngs modulus
% nu = Poissons ratio
% coord = coordinates at the element ends

c1=E*(1-nu)/(1+nu)/(1-2*nu);
c2=E*nu/(1+nu)/(1-2*nu);
c3=0.5*E/(1+nu);
C=[c1 c2 c2 0 0 0;
   c2 c1 c2 0 0 0;
   c2 c2 c1 0 0 0;
   0 0 0 c3 0 0;
   0 0 0  0 c3 0 ;
   0 0 0  0  0 c3];

[B,V]=bmat_tetra4(coord);
ke=B'*C*B*V;
