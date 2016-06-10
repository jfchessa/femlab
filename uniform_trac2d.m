function fext=uniform_trac2d(elem,node,trac,fext)

% function fext=uniform_trac2d(elem,node,trac,fext)
%
% Integrates a constant traction, trac, on a set of 1D boundary elements, 
% elem for a 2D problem.
%
% Works only for 2 node line elements at present
%


ne=size(elem,1);
nne=size(elem,2);
sdim=2;

for e=1:ne
  
  conne=elem(e,:);
  
  le = norm( node(elem(e,2),:) - node(elem(e,1),:) );
  
  fx = trac(1)*le;
  fy = trac(2)*le;
  
  fe = 0.5*[ fx; fy; fx; fy ];
  
  sctr(1:2:2*nne) = 2*conne-1;
  sctr(2:2:2*nne) = 2*conne;
  
  fext( sctr ) = fext( sctr ) + fe;
  
end
