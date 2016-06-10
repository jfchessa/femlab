% frame1.m
%
% A simple 2D frame problem solved using the finite element method
% in MATLAB.
%
% Written by Jack Chesssa, jfchessa@utep.edu
% Applied Finite Element Analysis
%
%

clear         % clear memory 

% ************************* INPUT SECTION *****************************

% define the node coordinate matrix
node = [  0.0  0.0;
         10.0  0.0;
         10.0  15.0;
          0.0  15.0 ];
  
% define the element connectivity matrix
connectivity = [ 2 3;
                 3 4;
                 4 1 ];

% define the element cross sectional area and Young's modulus
area = [ .1; .1; .1; .1 ];
young = 10e6*[ 1; 1; 1; 1]; 
Izz = 10e-3*[ 1; 1; 1; 1]; 
cpoints = [0.1; 0.1; 0.1; 0.1 ]; 

% define the global load vector
f=[ 0; 0; 0;   0; 0; 0;   5000; 0; 0;   0; 0; 0 ];

% define the fixed global dofs
ifix=[1 2 3 4 5 6];

% *********************** PROCESSING SECTION **************************

nn=size(node,1);            % number of nodes
ne=size(connectivity,1);    % number of elements
ndof=3*nn;                  % number of global dofs
K=sparse(ndof,ndof);        % define the global stiffness matrix

% assemble global stiffness matrix
for e=1:ne
    
    conn=connectivity(e,:);         % element connectivity
    coord=node(conn,:);             % element nodal coordinate matrix
    ke=kmat_beam2d(coord,Izz(e)*young(e),area(e)*young(e));  % get element stiffness matrix 
    
    sctr=get_scatter(conn,3); % get the global scatter vector 
    
    K(sctr,sctr)=K(sctr,sctr)+ke;   % scatter ke into K
    
end

% solve the system
[d,freac]=fesolve(K,f,ifix);


% ******************** POST-PROCESSING SECTION ************************

% write basic information about the model
fprintf('\n********************************************************');
fprintf('\n***         FINITE ELEMENT 2D FRAME PROGRAM          ***');
fprintf('\n********************************************************\n\n');
fprintf('number of nodes:      %4i\n',nn);
fprintf('number of elements:   %4i\n',ne);
fprintf('number of global dofs:%4i\n',ndof);

% write the nodal data
%x=reshape(d,2,3)'+node;     % the displaced nodal coordinates
fprintf('\nNodal Data\n');
fprintf('   NID    X-COORD     Y-COORD      X-DISP      Y-DISP        THETA  \n');
fprintf('--------------------------------------------------------------------\n');
for n=1:nn
    fprintf('%5i   %+8.3e   %+8.3e  %+8.3e  %+8.3e   %+8.3e  \n', n, node(n,1),...
        node(n,2), d(3*n-2), d(3*n-1), d(3*n) );    
end

% compute the element strains and stresses
fprintf('\nElement Data\n');
fprintf('   EID     CONNECTIVITY      STRAIN            STRESS\n');
fprintf('--------------------------------------------------------\n');
for e=1:ne
    
    conn=connectivity(e,:);
    n1=conn(1); n2=conn(2); 
    coord=node(conn,:);
    sctr=get_scatter(conn,3);
    
    B0=bmat_beam2d(coord,0);         % compute B-matrix at element midpoint
    B1=bmat_beam2d(coord,0);         % compute B-matrix at node 1
    B2=bmat_beam2d(coord,0);         % compute B-matrix at node 2
    y=cpoints(e);
    strainsb=[ B0*d(sctr)  B1*d(sctr) B2*d(sctr) ]*y;
    [val,indx]=max(abs(strainsb));
    strain(e)=strainsb(indx)+ bmat_truss2d(coord)*d([1 2 4 5]);        
    stress(e)=young(e)*strain(e);   % compute the element stress
    
    fprintf('%5i    %5i  %5i      %+10.5e     %+10.5e\n', e, n1, n2, ...
       strain(e), stress(e) );
    
    
end

% write the reaction force data
fprintf('\nReaction Force Data\n');
fprintf('  GDOF         FORCE\n');
fprintf('-------------------------\n');
for n=1:length(ifix)
    fprintf('%5i      %+10.5e\n', ifix(n), freac(n));
end


fprintf('\n*********************************************************\n\n');
