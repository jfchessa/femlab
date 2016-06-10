function F=projectInterfaceSpeed(node,conn,phi,zlsNode,speed)

% function F=projectInterfaceSpeed(node,conn,phi,zlsNode,speed)
%
% This funciton projects a interfacial speed defined on an iterface to
% the rest of the nodes in the computational domain.

numnode=size(node,1);
numelem=size(conn,1);

auglag=0;

% GET SYSTEM MATRIX A
A=sparse(numnode,numnode);
for e=1:numelem
  
  sctr=conn(e,:);
  [Q,W]=element_quadrature('tria3',1);
  
  % get tau
  h=sqrt(det([ [1;1;1]  node(sctr,:) ])/2);
  v=[ -1 1 0; -1 0 1 ]*phi(sctr);
  tau=h/abs(norm(v));
  
  for q=1:length(W)
    
    pt=Q(q,:);
    wt=W(q);
    N=shape_tria3(Q(q,:));
    dNdxi=dshape_tria3(Q(q,:));
    J=node(sctr,:)'*dNdxi;
    detJ=det(J); 
    dNdx=dNdxi*inv(J);
    
    phiPt=N'*phi(sctr);
    gradPhi=dNdx'*phi(sctr);
    
    % add gradPhi . gradF operator
    A(sctr,sctr)=A(sctr,sctr)+N*(gradPhi'*dNdx')*wt*detJ;
    % add SUPG operator
    A(sctr,sctr)=A(sctr,sctr)+(dNdx*gradPhi)*tau*(gradPhi'*dNdx')*wt*detJ;
    
  end % of quadrature loop
  
end % of element loop

% GET CONSTRIANT MATRIX G
nn=size(zlsNode,1);
G=sparse(nn,numnode);
l=zeros(nn,1);

% get elements zlsNodes are in
[es,xis] = tsearchn(node,conn(:,1:3),zlsNode);
edgeLookUp=[2 3;3 1;1 2];

% loop over zlsNodes
for n=1:nn
  
  % find edge cut by zls
  e=es(n);
  xi=xis(n,:);
  [val,ed]=min(xi);
  [ns]=edgeLookUp(ed,:);
  n1=conn(e,ns(1));
  n2=conn(e,ns(2));
  
  phi1=phi(n1);
  phi2=phi(n2);
  
  G(n,[n1 n2])=G(n,[n1 n2])+[phi2 -phi1]/(phi2-phi1);
  l(n)=speed(n);
  
end % loop over zlsNodes

% find nodes that are basically on the zls
zlsTol=1e-8;
ifix = find(abs(phi)<zlsTol);
ival = [];
if  ( length(ifix)>0 )
    
    disp('FOUND A NEAR NODE')
end

% solve system
beta=(10e3)*trace(A)/numnode;   

if ( auglag == 1 )
  A=[A+beta*G'*G    G';
     G              sparse(nn,nn)];   % augmented lagrangian
  l=[beta*G'*l; l];
else
  A=A+beta*G'*G;      % just a penalty method
  l=beta*G'*l;
end

F=fesolve(A,l,ifix,ival);
%F=A\l;
F=F(1:numnode);


