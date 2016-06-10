clear
fprintf('\n READING INPUT\n');

% ------------------------ DATA INPUT SECTION -------------------------
[node,nid]=readnodes('mode1.msh');       % read in mesh from
element=readelements('mode1.msh',102);   % gmsh file
nfix=readnodeset('mode1.msh',100);
nload=readnodeset('mode1.msh',101); 
                                        
node=node(:,1:2);     % renumber nod ids in conn, nfix and nload since
element=renumber(element,nid);  % gmsh does not always number them 
nfix=renumber(nfix,nid);        % consecutively
nload=renumber(nload,nid);


nn=size(node,1);  % number of nodes
ndof=2*nn;        % number of dofs
ne=size(element,1);  % number of elements

young = 10E6;
poisson = 0.33;
thk=1;

fext=zeros(ndof,1); 

ifix=[ 2*nfix-1 2*nfix 2*nload-1 ]; % fix nfix node in x and y
                                    % and nload node in x
fext(2*nload)=1000;   % point load at nload in the y-direction

fea2d_tria3