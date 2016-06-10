function F=vel2speed(node,conn,v,phi)

% function F=vel2speed(node,conn,v,phi)
%
% Takes a vector velocity field and a level set field on a finite element
% descretization and returns the scalar spped field F

numNode=size(node,1);
numElem=size(conn,1);

if ( size(v,1) == 2*numNode )
  v=[v(1:numNode) v(numNode+1:2*numNode)];
end

F=zeros(numNode,1);

A=sparse(numNode,numNode);
b=zeros(numNode,1);

[Q,W]=element_quadrature('tria3',2);
for e=1:numElem
  
  sctr=conn(e,:); 
   
  for q=1:size(W)
    
    N=shape_tria3(Q(q,:));
    dNdxi=dshape_tria3(Q(q,:));
    J=node(sctr,:)'*dNdxi;
    detJ=det(J);
    dNdx=dNdxi*inv(J);
    wt=W(q);
    
    vpt=N'*v(sctr,:);
    
    normal=dNdx'*phi(sctr);
    normal=normal/norm(normal);
    
    A(sctr,sctr)=A(sctr,sctr)+N*N'*detJ*wt;
    b(sctr)=b(sctr)+N*(vpt*normal)*detJ*wt;
    
  end % of quadrature loop
    
end % of element loop
  
F=A\b;  