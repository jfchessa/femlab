function [zlsNode,zlsConn]=lsmesh(node,conn,phi)

% function [zlsNode,zlsConn]=lsmesh(node,conn,phi)
%
% This function finds the zero level set mesh

numNode=size(node,1);   
numElem=size(conn,1);
eps=0;

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
denom=cutPhi(:,2)-cutPhi(:,3);
dz=find(denom==0);
denom(dz)=(cutPhi(dz,2)+cutPhi(dz,3))/100;
edge1Xi=(cutPhi(:,2)+cutPhi(:,3))./denom;

denom=cutPhi(:,3)-cutPhi(:,1);
dz=find(denom==0);
denom(dz)=(cutPhi(dz,3)+cutPhi(dz,1))/100;
edge2Xi=(cutPhi(:,3)+cutPhi(:,1))./denom;

denom=cutPhi(:,1)-cutPhi(:,2);
dz=find(denom==0);
denom(dz)=(cutPhi(dz,1)+cutPhi(dz,2))/100;
edge3Xi=(cutPhi(:,1)+cutPhi(:,2))./denom;

clear('dz');
clear('denom');

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

% eliminate redundant nodes
[zlsNode,inverseMap,newNumberMap]=unique(zlsNode,'rows');
zlsConn=newNumberMap(zlsConn);

% clf
% plot_mesh(node,conn,'T3','g')
% hold on
% plot_mesh(zlsNode,zlsConn,'L2','r')


