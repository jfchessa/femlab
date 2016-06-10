function phi=lsreinit(node,conn,phi)

% Function: lsreinit
%  
%       phi=lsreinit(node,conn,phi)
%
% This function takes reinitializes the level
% set function to be a signed distance function.
%
% by: Jack Chessa
%     Northwestern University
%     j-chessa@nwu.edu
%
%

N=size(node,1);   
E=size(conn,1);

if ( size(node,2) > 2 )
  node=node(:,1:2);
end

% *********************  MESH THE ZERO LEVEL SET ***************************
%%%%% find near nodes and create line mesh of zero level set and assign
%%%%% a value of F to each line element
[zlsNode,zlsConn]=lsmesh(node,conn,phi);

% get near nodes ( these are not to be changed )
cutElements=getCutElements(conn,phi);
nearNodes=unique(conn(cutElements,:));
farNodes=setdiff(1:N,nearNodes);

% *******************  CALCULATE PHI AT THE NODES **************************
% loop over the far nodes
for nf=1:length(farNodes)
  
  n=farNodes(nf);
  
  p=node(n,:);
  p1=zlsNode(zlsConn(:,1),:);
  p2=zlsNode(zlsConn(:,2),:);
  phi(n)=dist(p,p1,p2)*sign(phi(n));
  
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  