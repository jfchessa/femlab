node=[0 0;1 0;0 1];
sctr=[1 2 3];
phi=[-1;-1;1];
d=[1;1;1;1;1;1];

k1=1;
k2=10;

% discont quadrature
[W1,Q1]=discontT3quad(1,phi);


% brute force trapazoidal rule quadrature
ne=500;
[tri,X,Q2]=getSubElements2(node,ne);
W2=0.5/length(tri);



W=W1;Q=Q1;

K=zeros(6,6);
for q=1:length(Q)
  
  sctre=sctr+numnode;
  
  pt=Q(q,:);                             % quadrature point
  wt=W(q);                               % quadrature weight
  [N,dNdxi]=lagrange_basis('T3',pt);     % element shape functions
  
  % element Jacobian matrix w.r.t ref configuration
  J0=node(sctr,:)'*dNdxi;        
  invJ0=inv(J0);
  dNdx=dNdxi*invJ0;
  
  phiPt=N'*phis;
  signPhi=sign(phiPt);
  gradPhi=dNdx'*phis;
  
  psi=abs(phiPt)-abs(phis);
  gradPsi=sign(phiPt)*gradPhi;
  
  Ne=[N; (N.*filter).*psi];
  dNedx=[dNdx ; (dNdx.*[psi psi]+N*gradPsi').*[filter filter] ];
  
  if ( phiPt < 0 )
    k=k1;
  else
    k=k2;
  end
  
  K(sctre,sctre) = K(sctre,sctre) + dNedx'*dNedx*wt*det(J0);
  
end

K