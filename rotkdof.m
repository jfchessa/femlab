function K=rotkdof(K,R,gdof,ldof,nndof)

% function K=ROTKDOF(K,R,GDOF)
%
% Returns a the modified stiffness matrix and load vectors for locally 
% rotated (transformed) dofs.
%
%   K - the stiffness matrix or load vector
%   R - the rotation matrix
%   GDOF - a matrix of the global dofs to be rotated. Each row corresponds 
%       to the coupled dofs at a single node. This must have the same
%       number of coumns as the square dimension of R
%
%
% function K=ROTKDOF(K,R,NIDS,LDOF,NNDOF)
%
% Performs teh same function but when supplied node ids and local node dofs
% to be transformed. This assuems that the gdof spacing is as obtained by
% GET_SCATTER.
%
%   K - the stiffness matrix or load vector
%   R - the rotation matrix
%   LDOF - a vector for the local dofs to be transformed (one offset)
%   NNDOF -  the total number of dofs per node as used by GET_SCATTER.
%
%
% Part of FEMLAB
% Written by Jack Chessa, jfchessa@utep.edu  

n=size(K,1);
BigR=speye(n);

if ( nargin==3 )
    for i=1:size(gdof,1)
        BigR(gdof(i,:),gdof(i,:))=R; 
    end
else
    for i=1:length(gdof)
        ii=get_scatter( gdof(i), nndof, ldof );
        BigR(ii,ii)=R; 
    end
end

if ( size(K,2)==1 )
    K=BigR*K;
else
    K=BigR*K*BigR';
end

end