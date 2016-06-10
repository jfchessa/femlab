function nval=nodal_avg(eval,conn,node)

% function nval=nodal_avg(eval,conn,node)
%
% Computes nodal average values of element data
%

nn=size(node,1);
n=zeros(nn,1);

nval=zeros(nn,1);

if ( isstruct(conn) )
    for e=1:length(conn)
        
        sctr=conn(e).conn;
        nval(sctr)=nval(sctr)+eval(e);
        n(sctr)=n(sctr)+1;
        
    end
else
    
    for e=1:size(conn,1)
        
        sctr=conn(e,:);
        nval(sctr)=nval(sctr)+eval(e);
        n(sctr)=n(sctr)+1;
        
    end
end

nval=nval./n;
