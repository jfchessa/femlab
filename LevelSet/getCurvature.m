function kappa=getCurvature(node,element,gradPhi)

% function kappa=getCurvature(node,element,gradPhi)
%
% Computes curvature using a least squares projection

numNode=size(node,1);
numElem=size(element,1);

M=sparse(numNode,numNode);
f=zeros(numNode,1);

for e=1:numElem
  
  sctr=element(e,:);
  sctrx=sctr;
  sctry=sctr+numNode;
  
  [W,Q]=quadrature(2,'TRIANGULAR',2);
  for q=1:size(W,1)                      % quadrature loop
    
    pt=Q(q,:);                           % quadrature point
    wt=W(q);                             % quadrature weight
    [N,dNdxi]=lagrange_basis('T3',pt);   % element shape functions
    
    jac=node(sctr,:)'*dNdxi;  % element Jacobian matrix w.r.t ref configuration    
    dNdx=dNdxi*inv(jac);
    detj=det(jac);
    
    dm=N*N'*detj*wt;
    M(sctr,sctr)=M(sctr,sctr)+dm;
    
    divgrad=dNdx(:,1)'*gradPhi(sctr)+dNdx(:,2)'*gradPhi(sctr+numNode);
    a=divgrad/norm([ N'*gradPhi(sctr) N'*gradPhi(sctr+numNode) ]+eps);
    f(sctr)=f(sctr)+a*N*detj*wt;
    
  end
  
end 
kappa=M(1:numNode,1:numNode)\f;
