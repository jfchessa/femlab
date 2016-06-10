function [W,Q]=discontT6quad(order,phi1,phi2,level)

% function [W,Q]=discontT6quad(order,phi1,phi2,level)
%
% Returns an adaptive quadrature rule for quadratic triangular
% elements.  This quadrature subdivides the element so that it
% is refined around the discontinuity.

if ( nargin < 5 )
  level=4;
end
if ( nargin < 4 )
  phi2=phi1;
end

if ( max(phi1)*min(phi1)>=0  &  max(phi2)*min(phi2)>=0 )
  [W,Q]=quadrature(order,'TRIANGULAR',2);
  return;
end

% % DEBUG INPUT
% order=2;
% phi1=[-.4 -.4 .6 -.4 .1 .1]';
% phi1=[-.2 .7 .7 .1 .2 .1]';
% phi1=[.3 .3 .7 -.1 .4 .4]';
% phi2=phi1;
% level=4;

% ADAPTIVLY SUBDIVIDE ELEMENTS
subDivideNode=[ 0.00   0.00; 
                1.00   0.00; 
                0.00   1.00; 
                0.50   0.00; 
                0.50   0.50; 
                0.00   0.50; 
                0.25   0.00; 
                0.75   0.00; 
                0.75   0.25; 
                0.25   0.75; 
                0.00   0.75; 
                0.00   0.25; 
                0.50   0.25;  
                0.25   0.50; 
                0.25   0.25 ]; 
        
subDivideConn=[  1  4  6  7 15 12;
                 4  2  5  8  9 13;
                 6  5  3 14 10 11;
                 4  5  6 13 14 15  ];

subDivideN=zeros(15,6);
for n=1:15
  pt=subDivideNode(n,:);
  subDivideN(n,:)=lagrange_basis('T6',pt)';
end
   
subNode=subDivideNode;
subConn=subDivideConn;
subPhi1=subDivideN*phi1;
subPhi2=subDivideN*phi2;

for l=1:level 
  
  cutSubElem1=getCutElements(subConn,subPhi1);
  cutSubElem2=getCutElements(subConn,subPhi2);
  cutSubElem=unique([cutSubElem1;cutSubElem1]);
  
  % subdivide element
  nce=length(cutSubElem);
  for se=1:nce
    
    ne=size(subConn,1);
    nn=size(subNode,1);
    e=cutSubElem(se);
    sctr=subConn(e,:);
    
    % get new nodes
    nns=subDivideN*subNode(sctr,:);
    nns=nns(7:15,:);
    subNode=[subNode;nns];
    
    if ( se~=nce )
      cutSubElem(se+1:nce)=cutSubElem(se+1:nce)-1;
    end
    
    % add elements to connectivity
    subConn=subConn(setdiff(1:ne,e),:);
    subConn=[subConn;
             sctr(1) sctr(4) sctr(6) nn+1 nn+9 nn+6;
             sctr(4) sctr(2) sctr(5) nn+2 nn+3 nn+7;
             sctr(6) sctr(5) sctr(3) nn+8 nn+4 nn+5;
             sctr(4) sctr(5) sctr(6) nn+7 nn+8 nn+9 ];
             
    % get phi at new nodes
    newPhi1=subDivideN*subPhi1(sctr);
    newPhi2=subDivideN*subPhi2(sctr);
    subPhi1=[subPhi1;newPhi1(7:15)];
    subPhi2=[subPhi2;newPhi2(7:15)];
    
  end
  
%   clf
%   plot_field(subNode,subConn,'T6',subPhi1)
%   hold on
%   plot_mesh(subNode,subConn,'T6','k.-')
%   pause(0.1)
  
end

% MAP QUADRATURE POINTS IN SUB ELEMENTS
% loop over subtriangles to get quadrature points and weights
pt=1;
[w,q]=quadrature(order,'TRIANGULAR',2);
nq=length(w);
nse=size(subConn,1);
W=zeros(nse*nq,1);
Q=zeros(nse*nq,2);
for e=1:nse
  
  sctr=subConn(e,:);
  % transform quadrature points into the parent element
  coord=subNode(sctr,:);
  
  for n=1:length(w)
    [N,dNdxi]=lagrange_basis('T6',q(n,:));
    J=subNode(sctr,:)'*dNdxi;
    a=det(J)/2;
    Q(pt,:)=N'*coord;
    W(pt,1)=2*w(n)*a;
    pt=pt+1;
  end
 
end

% plot(Q(:,1),Q(:,2),'mx')