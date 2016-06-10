function [gradPhi,M,A]=getGradientField(node,element,phi,M,A)

% function [gradPhi,M,A]=getGradientField(node,element,phi,M,A)
%
% Computes a least squared projection for the level set gradient
%
% function gradPhi=getGradientField(node,element,phi)
% function [gradPhi,M,A]=getGradientField(node,element,phi)
% function gradPhi=getGradientField(node,element,phi,M,A)

numNode=size(node,1);
numElem=size(element,1);

if ( nargin == 3 ) % compute M and A
  
  M=sparse(2*numNode,2*numNode);
  A=sparse(2*numNode,numNode);
  
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
      M(sctrx,sctrx)=M(sctrx,sctrx)+dm;
      M(sctry,sctry)=M(sctry,sctry)+dm;
      A(sctrx,sctr)=A(sctrx,sctr)+N*dNdx(:,1)'*detj*wt;
      A(sctry,sctr)=A(sctry,sctr)+N*dNdx(:,2)'*detj*wt;
      
    end
    
  end 
end

f=A*phi;
gradPhi=M\f;
