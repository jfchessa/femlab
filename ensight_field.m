function stat=ensight_field(filename,data,etype)

% function stat=ensight_field(filename,data,etype)
%
% Write grid data to an ensight gold fromat.  If etype is given then
% it assumes per element data else per node data.
%
% etype Line2, Line3, Tria3, Tria6, Quad4, Quad8, 
%       Tetra4, Tetra10, Hexa8, Hexa20
%
% To print vector data have DATA be a matrix where the compnents are in row
% format.
%

fid=fopen(filename,'w');

if ( fid<0 ) 
  stat=1;
  return
else
  stat=0;
end

fprintf(fid,'%s\npart\n%10d\n',strtok(filename,'.'),1);
if ( nargin>2 )

  switch lower(etype)
    case 'line2'
      fprintf(fid,'bar2\n');
    case 'line3'
      fprintf(fid,'bar3\n');
    case 'tria3'
      fprintf(fid,'tria3\n');
    case 'tria6'
      fprintf(fid,'tria6\n');
    case 'quad4'
      fprintf(fid,'quad4\n');
    case 'quad8'
      fprintf(fid,'quad8\n');
    case 'tetra4'
      fprintf(fid,'tetra4\n');
    case 'tetra10'
      fprintf(fid,'tetra10\n');
    case 'hexa8'
      fprintf(fid,'hexa8\n');
    case 'hexa20'
      fprintf(fid,'hexa20\n');
    otherwise
      disp('Unknown etype in ENSIGHT_FIELD\n');
  end
  
else
  fprintf(fid,'coordinates\n');
end
fprintf(fid,'%12.5e\n', data);

fclose(fid);