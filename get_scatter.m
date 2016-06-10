function sctr=get_scatter( conn, nndof, ldof )

% function sctr=get_scatter( conn, nndof, ldof )
%
% Constructs a scatter vector, sctr, assuming a fixed number of dofs for
% each node.  The mapping of the local to global dofs, ldof is the global 
% dofs is as follows 
%
%   gdofi=(n-1)*nndof+ldofi
%
% where gdofi is a global dof, ldofi is a local dof and nndof is the number
% of dofs for each node.
%
%       conn - is the element connectivity
%       nndof - is the number of dofs for each node
%       ldof -  is a vector of the local dofs to map
%
%
% function sctr=get_scatter( conn, nndof )
%
% Works the same way but gives the mapping for all nndof for the element
% connectivity.
%
% Part of Femlab
% by Jack Chessa
% jfchessa@utep.edu
%

if ( nargin==2 )
    ldof=1:nndof;
end

lldof=length(ldof);
nn=length(conn);
lsctr=lldof*nn;

omat=(1:lldof)'*ones(1,nn);
emat=ones(lldof,1)*(nndof*(conn-1));

gmat=omat+emat;
sctr=reshape(gmat,1,lsctr);