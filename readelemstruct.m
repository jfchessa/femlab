function element=readelemstruct(meshFile,pids)


% function element=readelemstruct(filename)
%
%       Reads the element structures from a Gmsh ver 2.0 file
%
% function element=readelemstruct(filename,pids)
%
%       Reads the element structures from a Gmsh file with physical ids 
%       in pid.
%
% The element structure is of the following form
% element = struct('type',{'Tria3','Tria3','Quad4'},...
%   'pid',{1 1 1},...
%   'conn',{[1 2 3],[3 4 2],[1 2 5 6]});
% 
%  TYPE - element types: 
%       1  Line2    2-node line.
%       2  Tria3    3-node triangle.
%       3  Quad4    4-node quadrangle.
%       4  Tetra4   4-node tetrahedron.
%       5  Hexa8    8-node hexahedron.
%       6  Prism6   6-node prism.
%       7  Pyramid5 5-node pyramid.
%       8  Line3    3-node second order line  
%       9  Tria6    6-node second order triangle  
%      10  Quad9    9-node second order quadrangle  
%      11  Tetra10  10-node second order tetrahedron  
%      12  Hexa27   27-node second order hexahedron  
%      13  Prism18  18-node second order prism  
%      14  Pyramid14 14-node second order pyramid  
%      15  Point1   1-node point.
%      16  Quad8    8-node second order quadrangle  
%      17  Hexa20   20-node second order hexahedron  
%      18  Prism15  15-node second order prism  
%      19  Pyramid13 13-node second order pyramid  
%
% Written by Jack Chessa, jfchessa@utep.edu

% open the file
fid=fopen(meshFile,'r');
if ( fid<0 )
  disp(['Error could not open file ',meshFile]);
  node=[];
  nids=[];
  return
end

%* README: The 'msh' file format is the native output file format for
%  Gmsh. The file is divided in several sections (enclosed in $KEY and
%  $ENDKEY pairs). Two fields are important: $NOD/$ENDNOD defines the
%  nodes and $ELM/$ENDELM defines the elements.
%
%  The syntax is as follows:
%
%  $Nodes
%  number-of-nodes
%  node-number x-coord y-coord z-coord 
%  ...
%  $EndNodes
%
%  $Elements
%  number-of-elements
%  elm-number elm-type number-of-tags < tag > ... node-number-list
%  ...
%  $EndElements
%
%  All the syntactic variables stand for integers except x-coord,
%  y-coord and z-coord which stand for floating point values.  The
%  elm-type value defines the geometrical type for the element:
%
% number-of-tags
% gives the number of integer tags that follow for the n-th element. By default, the first tag is the number of the physical entity to which the element belongs; the second is the number of the elementary geometrical entity to which the element belongs; the third is the number of a mesh partition to which the element belongs. All tags must be postive integers, or zero. A zero tag is equivalent to no tag.
%
%  The elm-region value is the number of the physical entity to which
%  the element belongs. 
 
% read sections
while 1
  
  line=fgetl(fid);            % read line
  if ~isstr(line), break, end % check if EOF
  
  switch line                 % find seciton
    
  case '$Elements'          
    n=str2num(fgetl(fid));
    
    ne=0;
    for i=1:n
      temp=str2num(fgetl(fid));  % get element
      
      pid=temp(4);
      
      if ( nargin>1 )
        if ( nnz( ismember(pids,pid) ) == 0 )
          continue  % skip this element
        end
      end
      
      ne=ne+1;
      
      element(ne).eid=temp(1);
      etype=temp(2);
      ntags=temp(3);
      [estr,nn]=etypestr(etype);
      
      element(ne).pid=pid;
      element(ne).conn=temp(4+ntags:3+ntags+nn);
      element(ne).type=estr;
      
    end
    
  otherwise
    % skip line
  end
  
end

fclose(fid);



