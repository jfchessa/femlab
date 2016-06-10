function phi=flatten_levelset(node,conn,phi)

% function phi=flatten_levelset(node,conn,phi)
%
% This function flattens the level set 

numnode=size(node,1);
numelem=size(conn,1);

K=sparse(numnode,numnode);
fixedNodes=sparse(numnode,1);

for e=1:numelem
  
  sctr=conn(e,:);
  phis=phi(sctr);
  
  if ( max(phi(sctr))*min(phi(sctr)) < 0 )
    fixedNodes(sctr)=1;
  end
  
  [W,Q]=quadrature(1,'TRIANGULAR',2);      % define discontionus quadrature rule
  for q=1:length(W)                        % quadrature loop
    
    pt=Q(q,:);                             % quadrature point
    wt=W(q);                               % quadrature weight
    [N,dNdxi]=lagrange_basis('T3',pt);     % element shape functions
    
    % element Jacobian matrix w.r.t ref configuration
    J0=node(sctr,:)'*dNdxi;
    detJ0=det(J0);
    invJ0=inv(J0);
    dNdx=dNdxi*invJ0;
    
    K(sctr,sctr) =  K(sctr,sctr) + dNdx*dNdx'*wt*detJ0;
    
  end    
end

% solve
fixedNodes=find(fixedNodes);

f=-K(:,fixedNodes)*phi(fixedNodes);
f(fixedNodes)=phi(fixedNodes);
K(fixedNodes,:)=0;
K(:,fixedNodes)=0;
K(fixedNodes,fixedNodes)=speye(length(fixedNodes));

phi=K\f;


