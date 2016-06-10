function plot_mesh(X,connect,elem_type,se)

% function plot_mesh(X,connect,elem_type,linespec)
% 
% plots a nodal mesh and an associated connectivity.  X is
% teh nodal coordinates, connect is the connectivity, and
% elem_type is either 
%
%  'Line2', 'Line3', 'Tria3', 'Tria6', 'Quad4', 'Quad9', 'Hexa8' 
%
% depending on the element topology.
  
if ( nargin < 4 )
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
  
    switch lower(elem_type)
        case 'quad9'
            ord=[1,5,2,6,3,7,4,8,1];
            
        case 'tria3'
            ord=[1,2,3,1];
            
        case 'tria6'
            ord=[1,4,2,5,3,6,1];
            
        case 'quad4'
            ord=[1,2,3,4,1];
            
        case 'line2'
            ord=[1,2];
            
        case 'line3'
            ord=[1,3,2];
            
        case 'tetra4'
            ord=[1,2,4,1,3,4,2,3]; 
            
        case 'hexa8'
            ord=[1,5,6,2,3,7,8,4,1,2,3,4,8,5,6,7];
            
        otherwise
            disp(['Unknown element type ',elemt_type]);
    end
    
%    if ( strcmp(elem_type,'Quad9') )      % 9-node quad element
%       ord=[1,5,2,6,3,7,4,8,1];
%    elseif ( strcmp(elem_type,'Tria3') )  % 3-node triangle element
%       ord=[1,2,3,1];
%    elseif ( strcmp(elem_type,'Tria6') )  % 6-node triangle element
%       ord=[1,4,2,5,3,6,1];
%    elseif ( strcmp(elem_type,'Quad4') )  % 4-node quadrilateral element
%       ord=[1,2,3,4,1];
%    elseif ( strcmp(elem_type,'Line2') )  % 2-node line element
%       ord=[1,2];   
%    elseif ( strcmp(elem_type,'Line3') )  % 3-node line element
%       ord=[1,3,2];   
%    elseif ( strcmp(elem_type,'Tetra4') )  % 4-node tet element
%       ord=[1,2,4,1,3,4,2,3];   
%    elseif ( strcmp(elem_type,'Hexa8') )  % 8-node brick element
%       ord=[1,5,6,2,3,7,8,4,1,2,3,4,8,5,6,7];   
%    end
   
   for n=1:size(ord,2)
      xpt(n)=X(connect(e,ord(n)),1);
      ypt(n)=X(connect(e,ord(n)),2);      
      zpt(n)=X(connect(e,ord(n)),3);
   end
   plot3(xpt,ypt,zpt,se)
end

rotate3d on
axis equal
      
if ( ~holdState )
  hold off
end
