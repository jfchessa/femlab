function stat=ensight_fegeometry(filename,node,conn,etype)

% function stat=ensight_fegeometry(filename,node,conn)
%
% Writes finite element geometry to an Ensight Gold format file.  
%
%    filename - the name of the ensight file to be written (no extensions
%               are added)
%    node -  the node coordinate matrix
%    conn -  An array of element structure data (must contain .type and
%            .conn in the structure)
%
% function stat=ensight_fegeometry(filename,node,conn,etype)
%
% Writes finite element geometry to an Ensight Gold format file using a 
% standard element connectivity matrix. 
%
%    filename - the name of the ensight file to be written (no extensions
%               are added)
%    node -  the node coordinate matrix
%    conn -  A matrix of the element connectivities
%    etype - the element type.
%
% Element types supported are
%       etype ={ Line2, Line3, Tria3, Tria6, Quad4, Quad8, 
%                 Tetra4, Tetra10, Hexa8, Hexa27 }
%
%


fid=fopen(filename,'w');

if ( fid<0 )
    stat=1;
    return
else
    stat=0;
end

% write the nodes
if ( size(node,2) < 3 )
    node=[node,zeros(size(node,1),3-size(node,2))];
end

title=strtok(filename,'.');

fprintf(fid,'%s\n',filename);
fprintf(fid,'%s\n',title);
fprintf(fid,'%s\n','node id assign');
fprintf(fid,'%s\n','element id assign');
fprintf(fid,'%s\n','part');
fprintf(fid,'%10d\n', 1 );
fprintf(fid,'%s\n',title);

fprintf(fid,'%s\n','coordinates');
fprintf(fid,'%10d\n', size(node,1));
fprintf(fid,'%12.5e\n', node);

% write elements
if ( isstruct(conn) )

    es=1;
    ef=0;
    etype=conn(1).type;
    e=1;
    while ( es<=length(conn) )
        
        etype=conn(es).type;
        
        % count how many are of the same type
        while (  e<=length(conn) )
            if ( strcmp(conn(e).type,etype)==0 ), break;, end
            e=e+1;
        end
        ef=e-1;

        % go back and write that block
        [ename,connfmt,econn]=ensightconnfmt(etype);
        fprintf(fid,'%-80s\n',ename);
        fprintf(fid,'%10d\n', ef-es+1);
        for e=es:ef
            fprintf(fid,connfmt,conn(e).conn);
        end
        
        % increment
        es=ef+1;
        e=ef+2;

    end


else
%     switch (etype)
%         case 'Line2'
%             fprintf(fid,'%s\n','bar2');
%             fprintf(fid,'%10d\n', size(conn,1));
%             fprintf(fid,'%10d%10d\n',conn');
%         case 'Line3'
%             fprintf(fid,'%s\n','bar3');
%             fprintf(fid,'%10d\n', size(conn,1));
%             fprintf(fid,'%10d%10d%10d\n',conn');
%         case 'Tria3'
%             fprintf(fid,'%s\n','tria3');
%             fprintf(fid,'%10d\n', size(conn,1));
%             fprintf(fid,'%10d%10d%10d\n',conn');
%         case 'Tria6'
%             fprintf(fid,'%s\n','tria6');
%             fprintf(fid,'%10d\n', size(conn,1));
%             fprintf(fid,'%10d%10d%10d%10d%10d%10d\n',conn');
%         case 'Quad4'
%             fprintf(fid,'%s\n','quad4');
%             fprintf(fid,'%10d\n', size(conn,1));
%             fprintf(fid,'%10d%10d%10d%10d\n',conn');
%         case 'Quad9'
%             fprintf(fid,'%s\n','quad8');
%             fprintf(fid,'%10d\n', size(conn,1));
%             fprintf(fid,'%10d%10d%10d%10d%10d%10d%10d%10d\n',conn(:,1:8)');
%         case 'Tetra4'
%             fprintf(fid,'%s\n','tetra4');
%             fprintf(fid,'%10d\n', size(conn,1));
%             fprintf(fid,'%10d%10d%10d%10d\n',conn');
%         case 'Tetra10'
%             fprintf(fid,'%s\n','tetra10');
%             fprintf(fid,'%10d\n', size(conn,1));
%             fprintf(fid,'%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d\n',conn');
%         case 'Hexa8'
%             fprintf(fid,'%s\n','hexa8');
%             fprintf(fid,'%10d\n', size(conn,1));
%             fprintf(fid,'%10d%10d%10d%10d%10d%10d%10d%10d\n',conn');
%         case 'Hexa27'
%             fprintf(fid,'%80s\n','hexa20');
%             fprintf(fid,'%10d\n', size(conn,1));
%             fprintf(fid,'%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d\n',conn(:,1:20)');
%         otherwise
%             disp('unknown element type in ENSIGHT-FEGEOMETRY')
%     end

    [ename,connfmt,econn]=ensightconnfmt(etype);
    fprintf(fid,'%-80s\n',ename);
    fprintf(fid,'%10d\n', size(conn,1));
	fprintf(fid,connfmt,conn(:,econn)');
    
end

fclose(fid);

function [ename,connfmt,econn]=ensightconnfmt(etype)
    switch (lower(etype))
        case 'line2'
            ename='bar2';
            connfmt='%10d%10d\n';
            econn=1:2;
        case 'line3'
            ename='bar3';
            connfmt='%10d%10d%10d\n';
            econn=1:3;
        case 'tria3'
            ename='tria3';
            connfmt='%10d%10d%10d\n';
            econn=1:3;
        case 'tria6'
            ename='tria6';
            connfmt='%10d%10d%10d%10d%10d%10d\n';
            econn=1:6;
        case 'quad4'
            ename='quad4';
            connfmt='%10d%10d%10d%10d\n';
            econn=1:4;
        case 'quad9'
            ename='quad8';
            connfmt='%10d%10d%10d%10d%10d%10d%10d%10d\n';
            econn=1:8;
        case 'tetra4'
            ename='tetra4';
            connfmt='%10d%10d%10d%10d\n';
            econn=1:4;
        case 'tetra10'
            ename='tetra10';
            connfmt='%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d\n';
            econn=1:10;
        case 'hexa8'
            ename='hexa8';
            connfmt='%10d%10d%10d%10d%10d%10d%10d%10d\n';
            econn=1:8;
        case 'hexa27'
            ename='hexa20';
            connfmt='%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d\n';
            econn=1:20;
        otherwise
            disp('unknown element type in ENSIGHT-FEGEOMETRY');
            ename='';
            connfmt='\n';
            econn=[];
    end