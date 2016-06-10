function enriched=getEnrichedNodes(node,conn,phi)

% function enrichedNodes=getEnrichedNodes(node,conn,phi)

test=phi(conn');
test=max(test).*min(test);
cutElements=find(test < 0);

% find enriched nodes
enrichedNodes=unique(conn(cutElements,:));

enriched=sparse(size(node,1),1);
enriched(enrichedNodes)=1;