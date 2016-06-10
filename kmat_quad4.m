function ke = kmat_quad4(coord, E, nu, thk)

%  ke = kmat_quad4(coord, E, nu, thk)
% Generates equations for a plane stress 4 node quadralateral element, ke 
% E = Youngs modulus
% nu = Poissons ratio
% thk = element thickness
% coord = coordinates at the element ends

c1=E/(1-nu^2);  % plane stress material stiffness matrix
c2=nu*c1;
c3=0.5*(1-nu)*c1;
C=[c1 c2 0;c2 c1 0; 0 0 c3];

qpt=sqrt(1/3)*[ -1 -1;
                 1 -1;
                 1  1;
                 -1 1];
ke=zeros(8,8);
for q=1:4
	xi=qpt(q,:);
	[B,jac]=bmat_quad4(coord,xi);
	ke = ke + B'*C*B*jac*thk;
end
