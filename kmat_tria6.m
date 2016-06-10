function ke = kmat_tria6(coord, E, nu, thk)

%  ke = kmat_tria6(coord, E, nu, thk)
% Generates equations for a plane stress 6 node triangle element, ke 
% E = Youngs modulus
% nu = Poissons ratio
% thk = element thickness
% coord = coordinates at the element ends

c1=E/(1-nu^2);  % plane stress material stiffness matrix
c2=nu*c1;
c3=0.5*(1-nu)*c1;
C=[c1 c2 0;c2 c1 0; 0 0 c3];

qpt = [ 0.1666666666667, 0.1666666666667 ;
    0.6666666666667, 0.1666666666667 ;
    0.1666666666667, 0.6666666666667 ];

qwt = [0.3333333333333, 0.3333333333333, 0.3333333333333];

ke=zeros(12,12);
for q=1:3
  xi=qpt(q,:);
  [B,jac]=bmat_tria6(coord,xi);
  ke = ke + B'*C*B*jac*thk*qwt(q);
end
