% truss1.m
%
% A simple 2D truss problem solved using the finite element method
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
          0.0 15.0 ];
  
% define the element connectivity matrix
connectivity = [ 1 2;
                 2 3;
                 3 1 ];

% define the element cross sectional area and Young's modulus
area = [ .1; .15; .1 ];
young = 10e6*[ 1; 1; 1]; 

% define the global load vector
f=[ 0; 0; 0; -5000; 0; 0 ];

% define the fixed global dofs
ifix=[1 5 6];

% *********************** PROCESSING SECTION **************************

nn=size(node,1);            % number of nodes
ne=size(connectivity,1);    % number of elements
ndof=2*nn;                  % number of global dofs
K=sparse(ndof,ndof);        % define the global stiffness matrix

% assemble global stiffness matrix
for e=1:ne
    
    conn=connectivity(e,:);         % element connectivity
    coord=node(conn,:);             % element nodal coordinate matrix
    ke=kmat_truss2d(coord,area(e)*young(e));  % get element stiffness matrix 
    
%    sctr=get_scatter(conn,2,[1 2]); % get the global scatter vector 
    sctr=[ 2*conn(1)-1, 2*conn(1), 2*conn(2)-1, 2*conn(2) ]; % get the global scatter vector [u1 v1 u2 v2]
    K(sctr,sctr)=K(sctr,sctr)+ke;   % scatter ke into K
    
end

% solve the system
[d,freac]=fesolve(K,f,ifix);


% ******************** POST-PROCESSING SECTION ************************

% write basic information about the model
fprintf('\n********************************************************');
fprintf('\n***         FINITE ELEMENT 2D TRUSS PROGRAM          ***');
fprintf('\n********************************************************\n\n');
fprintf('number of nodes:      %4i\n',nn);
fprintf('number of elements:   %4i\n',ne);
fprintf('number of global dofs:%4i\n',ndof);

% write the nodal data
x=reshape(d,2,3)'+node;     % the displaced nodal coordinates
fprintf('\nNodal Data\n');
fprintf('   NID    X-COORD     Y-COORD      X-DISP      Y-DISP  \n');
fprintf('--------------------------------------------------------\n');
for n=1:nn
    fprintf('%5i   %+8.3e   %+8.3e  %+8.3e  %+8.3e  \n',n,node(2*n-1),...
        node(2*n),d(2*n-1),d(2*n));    
end

% compute the element strains and stresses
fprintf('\nElement Data\n');
fprintf('   EID     CONNECTIVITY      STRAIN            STRESS\n');
fprintf('--------------------------------------------------------\n');
for e=1:ne
    
    conn=connectivity(e,:);
    n1=conn(1); n2=conn(2); 
   
    le=abs(node(n2,:)-node(n1,:));  % undeformed length of the element
    lp=abs(x(n2,:)-x(n1,:));        % and the deformed length
    strain(e)=(lp-le)/le;           % compute the element strain
    stress(e)=young(e)*strain(e);   % compute the element stress
    
    fprintf('%5i    %5i  %5i      %+10.5e     %+10.5e\n',e,n1,n2,...
        strain(e),stress(e));
    
    
end

% write the reaction force data
fprintf('\nReaction Force Data\n');
fprintf('  GDOF         FORCE\n');
fprintf('-------------------------\n');
for n=1:length(ifix)
    fprintf('%5i      %+10.5e\n', ifix(n), freac(n));
end

% Plot the deformed truss
figure(1)
clf
hold on
title('DEFORMED PLOT OF TRUSS WITH ELEMENT STRESSES');
axis equal
nmappts=20;
cmap=jet(nmappts);
smin=min(stress); smax=max(stress);
srange=linspace(smin,smax,nmappts)';
for e=1:ne
    conn=connectivity(e,:);
    plot(node(conn,1),node(conn,2),'k--')
    ecolor=[ interp1(srange,cmap(:,1),stress(e)), ...
       interp1(srange,cmap(:,2),stress(e)), ...
       interp1(srange,cmap(:,3),stress(e))]; 
    plot(x(conn,1),x(conn,2),'LineWidth',2,'Color',ecolor)
end
caxis([smin,smax]);
colorbar

fprintf('\n*********************************************************\n\n');
