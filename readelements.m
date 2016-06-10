function [element,eid,type]=readelements(meshFile,pids)


% function [element,eid,type]=readelements(filename)
%
%       Reads the element connectivity matrix from a Gmsh ver 2.0 ASCII file
%
% function [element,eid,type]=readgmsh(filename,pids)
%
%       Reads the element connectivity matrix from a Gmsh file with 
%       physical ids in pid.
%
%
% This is part ot FEMLAB
% Written by Jack Chessa, jfchessa@utep.edu
%


% open the file
fid=fopen(meshFile,'r');
if ( fid<0 )
  disp(['Error could not open file ',meshFile]);
  element=[];
  eid=[];
  type=[];
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
      
      eid(ne)=temp(1);
      etype=temp(2);
      [estr,nn]=etypestr(etype);
      type{ne}=estr;
      ntags=temp(3);
      element(ne,1:nn)=temp(4+ntags:3+ntags+nn);
      
    end
    
  otherwise
    % skip line
  end
  
end

fclose(fid);

