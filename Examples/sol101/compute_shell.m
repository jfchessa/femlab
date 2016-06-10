function K=compute_shell( nodeCoord, dofmap, element, prop, mat, K )

etype='none';
sdim=3;
nndof=6;
for e=1:length( element )
    
    econn=element(e).conn;
    
    gcoord=nodeCoord(econn,:);
    thk=prop(element(e).pid).thk;
    mid=prop(element(e).pid).mid;
    young=mat(mid).young;
    nu=mat(mid).poisson;
    
    if ( etype!=element(e).type )
        nne=length(econn);
        ke=zeros(nne*nndof);
    else
        ke=0;
    end
    esctr=reshape( dofmap(econn,1:6), 1, nne*nndof );
    
    switch lower(etype)
        
        %case 'tria3'
        
        %case 'tria6'
        
        case 'quad4'
            ke=kmat_mindlin_quad4( gcoord, thk, young, nu );
            
        %case 'quad8'
            
        otherwise
            disp([etype,' element not supported in compute_shells']);
            
    end
    
    K(esctr,esctr)=K(esctr,esctr)+ke;
    
end
