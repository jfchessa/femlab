function status = write_abaqus(filename,node,conn,nset)
%  write_abaqus(filename,node,conn,pid,type)
%
%  Writes an Abaqus *.inp input file
% 
%  filename - 
%  node - 
%  conn -
%  nset - nodeset (in cell format if multiple sets)
%

%if ( length(pid)<size(conn,1) )
%    pid=pid(1)*ones(1,size(conn,1));
%end

nn=size(node,1);
while ( size(node,2)<3 )
    node = [node zeros(nn,1)];
end

%open the file
fid=fopen(filename,'w');

if ( fid<0 )
    stat=1;
    return
else
    stat=0;
end

% write the nodes
fprintf(fid,'**-------------------------------------------------------\n');
fprintf(fid,'**\n');
fprintf(fid,'**    N o d e     D e f i n i t i o n\n');
fprintf(fid,'**\n');
fprintf(fid,'**-------------------------------------------------------\n');
fprintf(fid,'**\n*Node\n');
for i=1:size(node,1)
    fprintf(fid,'%d, %f,%f,%f\n', i, node(i,1), node(i,2), node(i,3) );
end

if ( nargin<3 )
    return
end

% write the elements
fprintf(fid,'**\n**\n**-------------------------------------------------------\n');
fprintf(fid,'**\n');
fprintf(fid,'**    E l e m e n t     D e f i n i t i o n\n');
fprintf(fid,'**\n');
fprintf(fid,'**-------------------------------------------------------\n');
fprintf(fid,'**\n*Element, Type=sr4, Elset=name');
for i=1:size(conn,1)
    fprintf(fid,'\n%d', i);
    for j=1:size(conn,2)
        fprintf(fid,', %d', conn(i,j) );
    end
end


if ( nargin<4 )
    return
end

% write the node set(s)
fprintf(fid,'\n**\n**\n**-------------------------------------------------------\n');
fprintf(fid,'**\n');
fprintf(fid,'**    N o d e    S e t   D e f i n i t i o n s\n');
fprintf(fid,'**\n');
fprintf(fid,'**-------------------------------------------------------\n');
if ( iscell(nset) )
    for n=1:length(nset)
        fprintf(fid,'**\n*Nset, Nset=%s\n',['nset',num2str(n)]);
        r=0;
        for i=1:length(nset{n})
            fprintf(fid,'%d', nset{n}(i));
            r=r+1;
            if ( r==10 || i==length(nset{n}) )
                r=0;
                fprintf(fid,'\n');
            else
                fprintf(fid,', ');
            end
        end
    end
else
end
%*NSET, NSET=GEN1, GENERATE

fprintf(fid,'**\n**\n');

end

