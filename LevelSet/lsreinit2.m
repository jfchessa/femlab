function phi=lsreinitT3(node,conn,phi,eps,nsteps)

% Function: lsreinitT3
%  
%       phi=lsreinitT3(node,conn,phi,eps,nsteps)
%
% This function takes reinitializes the level
% set function to be a signed distance function.
%
% by: Jack Chessa
%     Northwestern University
%     j-chessa@nwu.edu
%
%

numNode=size(node,1);   
numElem=size(conn,1);

% get element size factors
sizeFact=getElementSizeFactor(node,conn,'T3');
Cr=0.9;

if ( nargin < 5)
  nsteps=5;
end
if ( nargin < 4)
  eps=2*min(sizeFact);
end

if ( size(node,2) > 2 )
  node=node(:,1:2);
end

% get elements not intersected by the zls
cutElem=getCutElements(conn,phi);
unCutElem=setdiff(1:numElem,cutElem);
integrateElem=1:numElem;

% vector of dofs that are to be updated
freeDof=setdiff( [1:numNode]', unique(conn(cutElem,:)) );
%freeDof=[1:numNode]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET LUMPED MASS MATRIX
M=sparse(numNode,numNode);
[W,Q]=quadrature( 4, 'TRIANGULAR' , 2 );
for e=integrateElem

  sctr=conn(e,:);
  phis=phi(sctr);
  
  for q=1:length(W)
    
    pt=Q(q,:);
    wt=W(q);
    [N,dNdxi]=lagrange_basis('T3',pt);
    J=node(sctr,:)'*dNdxi;
    detJ=det(J); 
    
    M(sctr,sctr)=M(sctr,sctr)+N*N'*detJ*wt; %sum(N*N'*detJ*wt)';
    
  end % of quadrature loop
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:nsteps
 
  dtau=10e10;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % GET FORCE VECTOR
  f=zeros(numNode,1);
  fsupg=zeros(numNode,1);
  [W,Q]=quadrature( 3, 'TRIANGULAR', 2 );
  for e=integrateElem
    
    sctr=conn(e,:);
    phis=phi(sctr);
    
    for q=1:length(W)
      
      pt=Q(q,:);
      wt=W(q);
      [N,dNdxi]=lagrange_basis('T3',pt);
      J=node(sctr,:)'*dNdxi;
      detJ=det(J); 
      dNdx=dNdxi*inv(J);
      
      phiPt=N'*phi(sctr);
      gradPhi=dNdx'*phi(sctr);
      
      sf=2*heavisider(phiPt/eps,2,2)-1;
      vPt=sf*gradPhi;
      %vPt=sf*gradPhi/(norm(gradPhi)+eps);
      
      % check for critical time step for next calc
      delTauPt=Cr*sizeFact(e)/(2*norm(vPt)+eps);
      if ( delTauPt < dtau )
        dtau=delTauPt;
      end
      
      % add fint
      f(sctr)=f(sctr)-N*(vPt'*gradPhi)*wt*detJ;
      
      % add fext
      f(sctr)=f(sctr)+N*sf*wt*detJ;
      
      % add SUPG operator
      fsupg(sctr)=fsupg(sctr)+(dNdx*vPt)*(vPt'*gradPhi-sf)*wt*detJ;
      
    end % of quadrature loop
  end
  
  % UPDATE PHI
  f=f-(dtau/2)*fsupg;
  dphi=sparse(numNode,1);
  dphi(freeDof)=dtau*(M(freeDof,freeDof)\f(freeDof));
  phi(freeDof)=phi(freeDof)+dphi(freeDof);
  
%   clf
%   plot_mesh(node,conn,'T6')
%   hold on
%   contura(node(:,1),node(:,2),phi,-10:0.5:10)
%   title(['REINITALIZATION STEP n=',num2str(n)])
%   pause(0.1)
 
end