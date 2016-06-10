function cutElements=getCutElements(conn,phi)

% function cutElemets=getCutElements(conn,phi)
%
% Finds the elements cut by the zero level set
%

test=phi(conn)';
test=max(test).*min(test);
cutElements=find( test < 0 )';