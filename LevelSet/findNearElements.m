function [nearNode,nearConn,nodeMap]=findNearELements(node,conn,phi,bandWidth)

% function [nearNode,nearConn,nodeMap]=findNearELements(node,conn,phi,bandWidth)
%
% Finds elements within a narrow band about the zero level set.

elemPhi=phi(conn);
minAbsElemPhi=min(abs(elemPhi'))';
nearElements=find(minAbsElemPhi<=bandWidth);

nearConn=conn(nearElements,:);
nodeMap=unique(nearConn);
nearNode=node(nodeMap,:);
inverseMap=sparse(size(node,1),1);
inverseMap(nodeMap)=[1:length(nearNode)]';
nearConn=inverseMap(nearConn);

