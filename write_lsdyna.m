function write_lsdyna(filename,node,conn,pid,type)

%  write_lsdyna(filename,node,conn,pid,type)
%
%  Writes an ls-dyna input file
%
%
%  type = SOLID, SHELL, BEAM, DISCRETE ..

if ( length(pid)<size(conn,1) )
    pid=pid(1)*ones(1,size(conn,1));
end

%open the file
fid=fopen(filename,'w');

if ( fid<0 )
    stat=1;
    return
else
    stat=0;
end


fprintf(fid,'$---1---->----2---->----3---->----4---->----5---->----6---->----7---->----8---->\n$\n');
fprintf(fid,'$\n$ LS-DYNA Input File\n$\n$   Model Description - \n$\n$\n$   Units:\n$\n');
fprintf(fid,'$---1---->----2---->----3---->----4---->----5---->----6---->----7---->----8---->\n$\n');

% write the nodes
fprintf(fid,'*NODE\n');
fprintf(fid,'$    NID         X-COORD         Y-COORD         Z-COORD      TC      RC\n');
for i=1:size(node,1)
    fprintf(fid,'%8d%16f%16f%16f\n', i, node(i,1), node(i,2), node(i,3) );
end

% write the elements
fprintf(fid,'$\n*ELEMENT_%s\n', upper(type));
if ( upper(type) == 'SOLID' )
    
    fprintf(fid,'$      EID       PID\n');
    fprintf(fid,'$       N1        N2        N3        N4        N5        N6        N7        N8');
    for i=1:size(conn,1)
        fprintf(fid,'\n%10d%10d\n', i, pid(i));
        for j=1:size(conn,2)
            fprintf(fid,'%10d', conn(i,j) );
        end
    end
    
else
    
    fprintf(fid,'$      EID       PID        N1        N2        N3        N4        N5        N6');
    for i=1:size(conn,1)
        fprintf(fid,'\n%10d%10d', i, pid(i));
        for j=1:size(conn,2)
            fprintf(fid,'%10d', conn(i,j) );
        end
    end
    
end
fprintf(fid,'\n$\n');
fprintf(fid,'$---1---->----2---->----3---->----4---->----5---->----6---->----7---->----8---->\n$\n');
