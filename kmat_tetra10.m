function ke = kmat_tetra10(coord, E, nu)

%  ke = kmat_tetra10(coord, E, nu)
% Generates equations for a ten node tetrahedral element, ke 
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
 
 qpt = [ 0.58541020  0.13819660  0.13819660;
   0.13819660  0.58541020  0.13819660;
   0.13819660  0.13819660  0.58541020;
   0.13819660  0.13819660  0.13819660];
 qwt = [1; 1; 1; 1]/4;

ke=zeros(30,30);
for q=1:length(qwt)
  xi=qpt(q,:);
  [B,jac]=bmat_tetra10(coord,xi);
  ke = ke + B'*C*B*jac*qwt(q);
end
