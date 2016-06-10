function plot_fefield(X,connect,d,elem_type,se)

% function plot_fefield(X,connect,d,elem_type,linespec)
% 
% plots a nodal mesh and an associated connectivity.  X is
% teh nodal coordinates, connect is the connectivity, and
% elem_type is either 
%
%  'Line2', 'Line3', 'Tria3', 'Tria6', 'Quad4', 'Quad9', 'Hexa8' 
%
% depending on the element topology.
  
if ( nargin < 5 )
   se='k-';
end

holdState=ishold;
hold on

% fill X if needed
if (size(X,2) < 3)
   for c=size(X,2)+1:3
      X(:,c)=[zeros(size(X,1),1)];
   end
end

for e=1:size(connect,1)
  
   if ( strcmp(lower(elem_type),'quad9') )      % 9-node quad element
      ord=[1,5,2,6,3,7,4,8,1];
   elseif ( strcmp(lower(elem_type),'tria3') )  % 3-node triangle element
      ord=[1,2,3,1];
   elseif ( strcmp(lower(elem_type),'tria6') )  % 6-node triangle element
      ord=[1,4,2,5,3,6,1];
   elseif ( strcmp(lower(elem_type),'quad4') )  % 4-node quadrilateral element
      ord=[1,2,3,4,1];
   elseif ( strcmp(lower(elem_type),'line2') )  % 2-node line element
      ord=[1,2];   
   elseif ( strcmp(lower(elem_type),'line3') )  % 3-node line element
      ord=[1,3,2];   
   elseif ( strcmp(lower(elem_type),'tetra4') )  % 4-node tet element
      ord=[1,2,4,1,3,4,2,3];   
   elseif ( strcmp(lower(elem_type),'hexa8') )  % 8-node brick element
      ord=[1,5,6,2,3,7,8,4,1,2,3,4,8,5,6,7];   
   end
   
   for n=1:size(ord,2)
      xpt(n)=X(connect(e,ord(n)),1);
      ypt(n)=X(connect(e,ord(n)),2);      
      zpt(n)=X(connect(e,ord(n)),3);
      cpt(n)=d(connect(e,ord(n)));
      
   end
   
   fill3(xpt,ypt,zpt,cpt,'edgecolor','none');
   
end

if ( ~strcmp( lower(se),'none' ) )
    plot_mesh(X,connect,elem_type,se)
end

rotate3d on
axis equal
      
if ( ~holdState )
  hold off
end
