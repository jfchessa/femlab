function kappa=getCurvature(node,element,bndry,phi)

% returns the least squared projection of the nodal curvatures
%

numNode=size(node,1);
numElem=size(element,1);

M=sparse(numNode,numNode);
f=zeros(numNode,1);

for e=1:numElem
  
  sctr=element(e,:);
  phis=phi(sctr);
  
  [W,Q]=quadrature(2,'TRIANGULAR',2);
  for q=1:size(W,1)                      % quadrature loop
    
    pt=Q(q,:);                           % quadrature point
    wt=W(q);                             % quadrature weight
    [N,dNdxi]=lagrange_basis('T3',pt);   % element shape functions
    
    jac=node(sctr,:)'*dNdxi;  % element Jacobian matrix w.r.t ref configuration    
    dNdx=dNdxi*inv(jac);
    detj=det(jac);
    
    gradPhiPt=dNdx'*phis;
    
    M(sctr,sctr)=M(sctr,sctr)+N*N'*detj*wt;
    %f(sctr)=f(sctr)-dNdx*gradPhiPt/norm(gradPhiPt)*detj*wt;
    f(sctr)=f(sctr)-dNdx*gradPhiPt*detj*wt;
    
  end
  
end 

% find T3 elements which correspond to bndry
xPt=node(:,1);
yPt=node(:,2);
testPt=[mean(xPt(bndry'))' mean(yPt(bndry'))'];
[es,xiT3s]=tsearchn(node,element,testPt);
fint=f;
% Boundary loop
for e=1:size(bndry,1)
  
  sctr=bndry(e,:);
  sctrT3=element(es(e),:);
  phis=phi(sctrT3);
  
  pt=0;                                    % quadrature point
  ptT3=xiT3s(e,1:2);
  wt=2;                                    % quadrature weight
  [N,dNdxi]=lagrange_basis('L2',pt);       % element shape functions
  [NT3,dNT3dxi]=lagrange_basis('T3',ptT3); % element shape functions
  
  jac=node(sctr,:)'*dNdxi;        % element Jacobian matrix w.r.t ref configuration 
  jacT3=node(sctrT3,:)'*dNT3dxi;  % element Jacobian matrix w.r.t ref configuration    
  dNT3dx=dNT3dxi*inv(jacT3);
  detj=norm(jac);
  
  normal=[jac(2) -jac(1)];
  normal=normal/norm(normal);
  
%   norm(normal)
%   xstart=N'*node(sctr,:);
%   xstop=xstart+normal;
%   arrow(xstart,xstop);
%   hold on
  
  gradPhiPt=dNT3dx'*phis;
  dphidn=normal*gradPhiPt;
  
  f(sctr)=f(sctr)+N*dphidn/norm(gradPhiPt)*detj*wt;
  
  
end 

% Solve for curvature kappa
kappa=M\f;


fb=f-fint;
clf
% plot curvature
trisurf(element,node(:,1),node(:,2),0*node(:,2),kappa);
%trisurf(element,node(:,1),node(:,2),phi);
view(2)
shading interp
hold on
trimesh(element,node(:,1),node(:,2),'Color',[.2 .2 .2])
colorbar

% hold on
% xs=xPt(bndry);
% ys=yPt(bndry);
% arrow([xs(:,1) ys(:,1)],[xs(:,2) ys(:,2)] )
