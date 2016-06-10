function [element,node,nodeMap]=renumber(element,nid,node)

% function [element,node,nodeMap]=renumber(element,nid,node)
%
% function element=renumber(element,nodeMap)
%
% Renumbers the nodes so that they are one offset and continous.  Nodes
% that are not referenced in element are removed.
%
%

if ( nargin==2 ) % just remap using the nodemap nid

    if ( isstruct(element) )  	% renumber an element structure 
   		for e=1:length(element)
       		element(e).conn = nid( element(e).conn )';
   		end
	
	else 						% renumber an array
	
		if ( size(element,1)==1 )  % if it is a row array we ned to keep it so
			element=nid(element)';
		else
			element=nid(element);
		end
	
	end

else  % we need to compute the nodemap

   	invNidMap(nid)=1:length(nid);
    nodeMap=zeros(size(nid));
    
	if ( isstruct(element) )  	% setup node map with a element structure
   		for e=1:length(element)
       		nodeMap( element(e).conn, 1 ) = -1;
   		end
	else  % it is a connectivity array
		for e=1:size(element,1)
			for i=1:size(element,2)
   				nodeMap( element(e,i), 1 ) = -1;
			end
		end
	end

   	iis=find( nodeMap<0 );
   	nodeMap(iis)=1:length(iis);
   	node=node(invNidMap(iis),:);

   	for e=1:length(element)
        if ( isstruct(element) )
            element(e).conn = nodeMap( element(e).conn )';
        else
            element(e,:) = nodeMap( element(e,:) )';
        end
   	end
end
