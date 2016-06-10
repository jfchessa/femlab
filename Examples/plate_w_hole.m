clear
fprintf('\n READING INPUT\n');

% ------------------------ DATA INPUT SECTION -------------------------
[node,nid]=readnodes('plate_w_hole.msh');       % read in mesh from
element=readelements('plate_w_hole.msh',12);    % gmsh file
nfixx=readnodeset('plate_w_hole.msh',10);
nfixy=readnodeset('plate_w_hole.msh',9);
eload=readelements('plate_w_hole.msh',11);
                                        
node=node(:,1:2);     % renumber nod ids in conn, nfix and nload since
element=renumber(element,nid);  % gmsh does not always number them 
nfixx=renumber(nfixx,nid);      % consecutively
nfixy=renumber(nfixy,nid);     
eload=renumber(eload,nid);

nn=size(node,1);  % number of nodes
ndof=2*nn;        % number of dofs
ne=size(element,1);  % number of elements

young = 10E6;
poisson = 0.33;
thk=0.5;

fext=zeros(ndof,1); 
fext=uniform_trac2d(eload,node,[0 100],fext);

ifix=[ 2*nfixx-1; 2*nfixy ]; 

fea2d_tria3