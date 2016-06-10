function stat=ensight_blockgeometry(filename,node,nx,ny,nz)

% function stat=ensight_blockgeometry(filename,node,nx,ny,nz)
%
% Writes structured grid geometry to an Ensight Gold format file.  
%
%    filename - the name of the ensight file to be written (no extensions
%               are added)
%    node -  the node coordinate matrix
%    nx, ny, nz - the number of nodes in the x, y and z directions,
%                 respectively
%
%    in the case the nodes are given in i-j-k order.  In otherwords the
%    nth node in the node list corresponds to (k-1)*nx*ny+(j-1)*nx+i (loop
%    over the i nodes then the j then the k)
%
% function stat=ensight_blockgeometry(filename,X,Y,Z) 
% 
%    in this case the node coordinates are given in the matricies X,Y,Z 
%    as would be constructed in  MESHGRID  
%    


fid=fopen(filename,'w');

if ( fid<0 )
    stat=1;
    return
else
    stat=0;
end


% write to an Ensight Gold file

fprintf(fid,'Ensight Gold Geometry File\nOrdered grid geometry data\n');
fprintf(fid,'node id off\nelement id off\n');
fprintf(fid,'part\n1\nBlock grid\n');
fprintf(fid,'block\n');

if ( nargin == 4 )
    fprintf(fid,'%10d%10d%10d\n',size(node,1),size(node,2),size(node,3));
    
else
    fprintf(fid,'%10d%10d%10d\n',nx,ny,nz);
    fprintf(fid,'%12.5e\n',nodes);
end

fclose(fid);
