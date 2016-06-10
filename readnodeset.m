function nids=readnodeset(meshFile,pids)


% function nids=readnodeset(meshFile,pids)
%
%       Reads the nodes ids belonging to physical id pid from a Gmsh file
%
% Written by Jack Chessa, jfchessa@utep.edu
%

% element = struct('type',{'Tria3','Tria3','Quad4'},...
%   'pid',{1 1 1},...
%   'conn',{[1 2 3],[3 4 2],[1 2 5 6]});
% 
% node = [0 0 0;
%         1 0 0;
%         1 1 0;
%         0 1 0;
%         1 2 0;
%         0 2 0];
%       
% nids=[1 2 3 4 5 6];


% open the file
fid=fopen(meshFile,'r');
if ( fid<0 )
  disp(['Error could not open file ',meshFile]);
  nids=[];
  return
end


 
% read sections

conn=readelements(meshFile,pids);
nids=unique(conn);

fclose(fid);
