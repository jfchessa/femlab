function phi=lsupdate(node,conn,phi,v,delt)

% function phi=lsupdate(node,conn,phi,v,delt)
%
% This function performs an update of a level set field phi on a finite
% element descretization using a Taylor-Galerkin stabilized finite element
% scheme.  Node is the nodal coordinate matrix, conn is the element 
% connectivity matrix, F is the nodal speed function, delt is the 
% time step and elemType it the element type.

% INITIALIZE DATA STRUCTURES
numNode=size(node,1);
numElem=size(conn,1);
M=zeros(numNode,1);
fint=zeros(numNode,1);

% COMPUTE invM and fint
[Q,W]=element_quadrature('tria3',2);
for e=1:numElem

  sctr=conn(e,:); 
   
  for q=1:size(W)
    
    N=shape_tria3(Q(q,:));
    dNdxi=dshape_tria3(Q(q,:));
    J=node(sctr,:)'*dNdxi;
    detJ=det(J);
    dNdx=dNdxi*inv(J);
    
    gradPhi=phi(sctr)'*dNdx;
    
    vPt=N'*v(sctr,:);
    
    % advection part
    fint(sctr) = fint(sctr) - N*(vPt*gradPhi')*W(q)*detJ;  
    
    % stabilization part
    fint(sctr) = fint(sctr) - delt/2*(dNdx*vPt')*(vPt*gradPhi')*W(q)*detJ;
     
    % lumped mass
    M(sctr) = M(sctr) + sum(N*N'*W(q)*detJ)';
    
  end % of quadrature loop
  
end % of element loop

% UPDATE LEVEL SET
phi = phi + delt*fint./M;