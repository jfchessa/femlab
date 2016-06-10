function fext = compute_force( loadcase, dofmap, forcedata, fext )


   for lc=1:length(loadcase)
      
       lsid=loadcase(lc).loadid;
       
       gdofs = dofmap( loadcase(lc).nid, 1:6 );
       fval=loadcase(lc).value';
       ii=find(gdofs);
       nl=length(ii);
       fext( gdofs(ii), lc ) = reshape( loadcase(lc).value(ii), nl, 1 );
       
   end

end

