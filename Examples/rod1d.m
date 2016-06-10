%--------------------------------------------------------------------------
% rod1d.m
%
% A finite element code for 1D rod elements
%
% This command removes all variables from MATLAB's memory
clear

%--------------------------------------------------------------------------
% Input section
%--------------------------------------------------------------------------
% Now we define the following finite element mesh 
%   
%   x=0.0                          x=2.0
%   |>1======2=======3=======4=======5 ---> F=100
%       e=1     e=2     e=3     e=4
%
% In this problem we will use a four element mesh (five nodes).  First
% we define the nodal coordinates (x-coords) as a column vector called
% NODE as follows:
node=[ 0.0; 0.5; 1.0; 1.5; 2.0 ];

% Next we need to define the element connectivity CONN.   Each row defines
% the element connectivity.  Element 1 connects nodes 1 to 2, element 2 
% connects nodes 2 to 3 ... etc. 
conn=[ 1 2; 2 3; 3 4; 4 5 ];

% Finally, for this mesh all the elements have the same cross-sectional 
% area A and Young's modulus E.  We define these below
A=1.0;     % element cross sectional area
E=10e6;    % element Young's modulus

% Now we need to define the boundary conditions for the problem
% 
% first let's define the DOFs that are to be fixed in this problem.  For 
% this case this is just node 1 (or dof 1 since there is one dof/node)
ifix=[1];

% Next we define the externally applied nodal force vector F.  All external
% forces are zero except for node number 5 which has a positive force of 
% 100.  Note that this mush be a column vector.
f=[ 0.0; 0.0; 0.0; 0.0; 100.0 ];

%--------------------------------------------------------------------------
% This is the end of the input section.  There rest of the program solves
% the finite element problem.
%--------------------------------------------------------------------------

% Here we just define some parameters that are useful in the computation.
nn=size(node,1);     % the number of nodes in the mesh
ne=size(conn,1);     % the number of elements in the mesh
ndof=nn;             % the number of DOFs in the mesh    
K=sparse(ndof,ndof); % define the global stiffness matrix K

% Here is the meat of the program.  This is where we loop over the 
% individual elements, compute the element stiffness matrix and scatter it 
% into the global stiffness matrix K.
for e=1:ne   % loop over the elements
    
    sctr = conn(e,:);      % A vector of gdofs where the ke is scattered
	n1=conn(e,1);          % The first and second node ids in 
    n2=conn(e,2);          % the element connectivity
	x1=node(n1);           % The x-coordinates of the nodes in the element
    x2=node(n2);                           
	L=abs(x2-x1);                  % the length of the element
	ke=A*E/L*[1.0 -1.0;-1.0 1.0];  % the element stiffness matrix
    K(sctr,sctr)=K(sctr,sctr)+ke;  % and scatter ke into K
    
end % end of the element loop
    
% Now we solve the discrete finite element system of equations Kd=f for the
% nodal displacements d and the reaction forces freac using the function 
% FESOLVE. This is not a built in MATLAB function but a user defined one.
[d,freac]=fesolve(K,f,ifix);

%--------------------------------------------------------------------------
% Postprocessing section
%--------------------------------------------------------------------------
% While the node displacements are useful we also would like to know the 
% element stresses and strains.  We compute those now

% print out the nodal coordiantes and displacements
fprintf('\n------------------- R E S U L T S ---------------------\n')
fprintf('\nNID              X-COORD                  DISPLACEMENT\n')
fprintf('--------------------------------------------------------\n')
fprintf(' %d             %e                %e\n',[(1:nn); node'; d'])

% loop over the elements and compute the element strain and stress

fprintf('\n\nEID    NODE1 NODE2         STRAIN             STRESS  \n')
fprintf('--------------------------------------------------------\n')
for e=1:ne  

	n1=conn(e,1); n2=conn(e,2); % (This is the same as in the prior 
	x1=node(n1);  x2=node(n2);  %  element loop)
	L=abs(x2-x1);                
	d1=d(n1); d2=d(n2);         % the displacements of the nodes
	strain=(d2-d1)/L;  % compute the strain as change in length over length
    stress=E*strain;   % and the stress using Hooke's law - stress=E*strain   
    fprintf(' %d       %d     %d        %e        %e\n', e, n1, n2, strain, stress)
    
end % end of the element loop

