function [node,nids]=readnodes(meshFile)

% function [node,nids]=readnodes(filename)
%
%       Reads the nodes and node ids from a Gmsh file
%
% function [node]=readgmsh(filename)
%
%       Reads the nodes and elements from a Gmsh 2.0 ASCII file with 
%       physical ids in pid.  
%
%   filename - string with the filename of the gmsh file.  It must contain
%   the relative path as well as the extension, typically .msh
%
% Written by Jack Chessa, jfchessa@utep.edu
%


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
%  elm-type: 
%  
% 1  2-node line.
% 2  3-node triangle.
% 3  4-node quadrangle.
% 4  4-node tetrahedron.
% 5  8-node hexahedron.
% 6  6-node prism.
% 7  5-node pyramid.
% 8  3-node second order line (2 nodes associated with the vertices and 1 with the edge).
% 9  6-node second order triangle (3 nodes associated with the vertices and 3 with the edges).
% 10 9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face).
% 11 10-node second order tetrahedron (4 nodes associated with the vertices and 6 with the edges).
% 12 27-node second order hexahedron (8 nodes associated with the vertices, 12 with the edges, 6 with the faces and 1 with the volume).
% 13 18-node second order prism (6 nodes associated with the vertices, 9 with the edges and 3 with the quadrangular faces).
% 14 14-node second order pyramid (5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face).
% 15 1-node point.
% 16 8-node second order quadrangle (4 nodes associated with the vertices and 4 with the edges).
% 17 20-node second order hexahedron (8 nodes associated with the vertices and 12 with the edges).
% 18 15-node second order prism (6 nodes associated with the vertices and 9 with the edges).
% 19 13-node second order pyramid (5 nodes associated with the vertices and 8 with the edges).
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

    
  case '$Nodes'           
    numnode=str2num(fgetl(fid));
    node=zeros(numnode,3);
    nids=zeros(numnode,1);
    
    for i=1:numnode
      nodeline=str2num(fgetl(fid));
      node(i,:)=nodeline(2:4);
      nids(i)=nodeline(1);
    end

  otherwise
    % skip line
  end
  
end

fclose(fid);
