%function F=projectInterface(node,conn,phi,zlsNode,speed)

% function F=projectInterface(node,conn,phi,zlsNode,speed)
%
% This funciton projects a interfacial speed defined on an iterface to
% the rest of the nodes in the computational domain.
clear
clf

[node,connectivities,elemType]=msh2mlab('square.msh'); 
node=node(:,1:2);
[connectivities,node]=t3tot6(connectivities,node);
conn=connectivities{9};

numNode=size(node,1);
numElem=size(conn,1);
numStep=100;
delt=0.025;

x=node(:,1);
y=node(:,2);
r=2;
x0=10/2;
y0=1.25*r;

phi=r-sqrt((x-x0).^2+(y-y0).^2);

for n=1:5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numnode=size(node,1);
numelem=size(conn,1);

eps=.025;

% GET SYSTEM MATRIX A
A=sparse(numnode,numnode);
G=sparse(numnode,numnode);
l=sparse(numnode,1);
activeDof=sparse(1,numnode);
for e=1:numelem
  
  sctr=conn(e,:);
  phis=phi(sctr);
  if ( max(phis)*min(phis)<0 )
    activeDof(sctr)=1;
  end
  
  [W,Q]=discontT6quad( 5, phis, 5 );
  for q=1:length(W)
    
    pt=Q(q,:);
    wt=W(q);
    [N,dNdxi]=lagrange_basis('T6',pt);
    J=node(sctr,:)'*dNdxi;
    detJ=det(J); 
    dNdx=dNdxi*inv(J);
    
    phiPt=N'*phi(sctr);
    gradPhi=dNdx'*phi(sctr);
    
    hf=2*heavisider(phiPt,eps,2,1)-1;
    vPt=hf*gradPhi;
    tau=detJ/(vPt'*vPt);
    
    % add gradPhi . gradF operator
    A(sctr,sctr)=A(sctr,sctr)+N*(vPt'*dNdx')*wt*detJ;
    
    % add SUPG operator
    A(sctr,sctr)=A(sctr,sctr)+(dNdx*vPt)*tau*(vPt'*dNdx')*wt*detJ;
    
    % add interface constraint operator
    df=diracr(phiPt,eps,1,1);
    fPt=N'*node(sctr,2);
    G(sctr,sctr)=G(sctr,sctr)+df*N*N'*wt*detJ;
    l(sctr)=l(sctr)+df*N*fPt*wt*detJ;
    
  end % of quadrature loop
  
end % of element loop

% solve system
beta=10e3;
S=[A+beta*G'*G  G;
   G'        sparse(numnode,numnode)];
l=[beta*l;l]; 

activeDof=[1:numnode find(activeDof)+numnode];

F=sparse(2*numnode,1);
F(activeDof)=S(activeDof,activeDof)\l(activeDof);
F=F(1:numnode);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf
plot_field(node,conn,'T6',F)
hold on
plot_mesh(node,conn,'T6','g')
contura(node(:,1),node(:,2),phi,[0 0],'k-')
colorbar
pause(0.1)

phiOld=phi;
phi=lsupdateT6(node,conn,phi,F,.1); 
phi=lsreinitT6(node,conn,phi);
clf
plot_field(node,conn,'T6',phi)
colorbar
pause(0.1)

end

