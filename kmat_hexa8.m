function ke = kmat_hexa8(coord, E, nu)

%  ke = kmat_hexa8(coord, E, nu)
% Generates equations for a eight node hexahedral element, ke 
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
 
 qpt = [   0.57735026918963   0.57735026918963   0.57735026918963;
   0.57735026918963   0.57735026918963  -0.57735026918963;
   0.57735026918963  -0.57735026918963   0.57735026918963;
   0.57735026918963  -0.57735026918963  -0.57735026918963;
  -0.57735026918963   0.57735026918963   0.57735026918963;
  -0.57735026918963   0.57735026918963  -0.57735026918963;
  -0.57735026918963  -0.57735026918963   0.57735026918963;
  -0.57735026918963  -0.57735026918963  -0.57735026918963 ];

qwt=ones(8,1);

ke=zeros(24,24);
for q=1:length(qwt)
  xi=qpt(q,:);
  [B,jac]=bmat_hexa8(coord,xi);
  ke = ke + B'*C*B*jac;
end
