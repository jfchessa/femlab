function [zlsNode,zlsConn]=lsmesh(node,conn,phi)

% function [zlsNode,zlsConn]=lsmesh(node,conn,phi)
%
% This function finds the zero level set mesh

numNode=size(node,1);   
numElem=size(conn,1);
eps=1e-8;

%% find the elements cut by the level set and near nodes
cutElements=getCutElements(conn,phi);
numCutElem=length(cutElements);

cutConn=conn(cutElements,1:3);
xConn=node(:,1);
xConn=xConn(cutConn);
yConn=node(:,2);
yConn=yConn(cutConn);
cutPhi=phi(cutConn);

clear('cutConn');
clear('cutElements');

% find xi for intersection of level set with edge
edge1Xi=cutPhi(:,2)-cutPhi(:,3);
edge2Xi=cutPhi(:,3)-cutPhi(:,1);
edge3Xi=cutPhi(:,1)-cutPhi(:,2);

edge1Xi(find(edge1Xi==0))=eps;
edge2Xi(find(edge2Xi==0))=eps;
edge3Xi(find(edge3Xi==0))=eps;

edge1Xi=(cutPhi(:,2)+cutPhi(:,3))./edge1Xi;
edge2Xi=(cutPhi(:,2)+cutPhi(:,3))./edge2Xi;
edge3Xi=(cutPhi(:,2)+cutPhi(:,3))./edge3Xi;

% get x and y coords of edge intersection
edge1X=xConn(:,2).*(1-edge1Xi)/2+xConn(:,3).*(1+edge1Xi)/2;
edge1Y=yConn(:,2).*(1-edge1Xi)/2+yConn(:,3).*(1+edge1Xi)/2;

edge2X=xConn(:,3).*(1-edge2Xi)/2+xConn(:,1).*(1+edge2Xi)/2;
edge2Y=yConn(:,3).*(1-edge2Xi)/2+yConn(:,1).*(1+edge2Xi)/2;

edge3X=xConn(:,1).*(1-edge3Xi)/2+xConn(:,2).*(1+edge3Xi)/2;
edge3Y=yConn(:,1).*(1-edge3Xi)/2+yConn(:,2).*(1+edge3Xi)/2;

% filter out xis such that -1 <= xi < 1 ( edge is not intersected )
goodEdge1Elems=intersect(find(edge1Xi >= -1 ),find(edge1Xi < 1 ));
goodEdge2Elems=intersect(find(edge2Xi >= -1 ),find(edge2Xi < 1 ));
goodEdge3Elems=intersect(find(edge3Xi >= -1 ),find(edge3Xi < 1 ));

clear('edge1Xi');
clear('edge2Xi');
clear('edge3Xi');

% initial zlsNode
zlsNode=(-1e10)*ones(3*numCutElem,2);
zlsNode(3*goodEdge1Elems-2,1)=edge1X(goodEdge1Elems);
zlsNode(3*goodEdge1Elems-2,2)=edge1Y(goodEdge1Elems);
zlsNode(3*goodEdge2Elems-1,1)=edge2X(goodEdge2Elems);
zlsNode(3*goodEdge2Elems-1,2)=edge2Y(goodEdge2Elems);
zlsNode(3*goodEdge3Elems,1)=edge3X(goodEdge3Elems);
zlsNode(3*goodEdge3Elems,2)=edge3Y(goodEdge3Elems);

goodNodes=find(zlsNode(:,1)~=(-1e10));
zlsNode=zlsNode(goodNodes,:);

% initial zls connectivity
zlsConn=[1:2:2*numCutElem-1;2:2:2*numCutElem]';

% clf
% plot_mesh(node,conn,'T3','g')
% hold on
% plot_mesh(zlsNode,zlsConn,'L2','r')

% find elements with an edge coinsident with the zls
zphi1 = find(phi(conn(:,1))==0);
zphi2 = find(phi(conn(:,2))==0);
zphi3 = find(phi(conn(:,3))==0);

commEdge1 = intersect(zphi1,zphi2);
commEdge2 = intersect(zphi2,zphi3);
commEdge3 = intersect(zphi3,zphi1);

commConn = [conn(commEdge1,[1 2]);conn(commEdge2,[2 3]);conn(commEdge3,[3 1])];
commNode = node(commConn,:);
commConn = [1:size(commConn,1); (1:size(commConn,1))+size(commConn,1)]';

% add them to zlsNode and zlsConn
zlsConn=[zlsConn;commConn+size(zlsConn,1)];
zlsNode=[zlsNode;commNode];

% eliminate redundant nodes
n=1;
del=1e-10;
nn=size(zlsNode,1);
while n<=nn
  
  nnn=(n+1):nn;
  dist=[zlsNode(n,1)-zlsNode(nnn,1) zlsNode(n,2)-zlsNode(nnn,2)];
  norm=sqrt(dist(:,1).^2+dist(:,2).^2);
  
  dupn=find(norm<=del)+n;  % duplicate nodes
  
  nnmap=zeros(nn,1);
  nnmap(dupn)=n;
  nnmap( setdiff([1:nn]',dupn) ) = [ 1:(nn-length(dupn)) ]';
  
  zlsNode=zlsNode(setdiff([1:size(zlsNode,1)]',dupn),:);  
  zlsConn = nnmap(zlsConn);
  
  n=n+1;
  nn=size(zlsNode,1);
  
end

% clf
% plot_mesh(node,conn,'T3','k-')
% hold on
% plot_mesh(zlsNode,zlsConn,'L2','ro-')
% pause(0.5)


