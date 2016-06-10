function ke = kmat_quad8(coord, E, nu, thk)

%  ke = kmat_quad8(coord, E, nu, thk)
% Generates equations for a plane stress eight node quad element, ke 
% E = Youngs modulus
% nu = Poissons ratio
% thk = element thickness
% coord = coordinates at the element ends

c1=E/(1-nu^2);  % plane stress material stiffness matrix
c2=nu*c1;
c3=0.5*(1-nu)*c1;
C=[c1 c2 0;c2 c1 0; 0 0 c3];

qpt = [ -0.774596669241483, -0.774596669241483;
         0.0,               -0.774596669241483;
         0.774596669241483, -0.774596669241483;
        -0.774596669241483,  0.0;
         0.0,                0.0;
         0.774596669241483,  0.0;
        -0.774596669241483,  0.774596669241483;
         0.0,                0.774596669241483;
         0.774596669241483,  0.774596669241483 ];
         
qwt = [ 0.30864197530864 0.49382716049383 0.30864197530864 ...
        0.49382716049383 0.79012345679012 0.49382716049383 ...
        0.30864197530864 0.49382716049383 0.30864197530864 ];

ke=zeros(16,16);
for q=1:9
  xi=qpt(q,:);
  [B,jac]=bmat_tria6(coord,xi);
  ke = ke + B'*C*B*jac*thk*qwt(q);
end
