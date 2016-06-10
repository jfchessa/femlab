%----------------------------------------------------------------------
%
% fea2d.m
%
% A two dimensional finite element code 
%
% Written by Jack Chessa jfchessa@utep.edu
% for CE 5307 Theory of Finite Element Analysis
% September 6, 2007
%
%----------------------------------------------------------------------
clear

% ------------------------ DATA INPUT SECTION -------------------------

node=[0.0 0.0; 
      1.0 0.0; 
      2.0 0.0; 
      3.0 0.0; 
      0.5 1.0;
      1.5 1.0; 
      2.5 1.0];  % node coordinate matrix
    
conn=[1 2;
      2 3;
      3 4;
      5 6;
      6 7;
      1 5;
      5 2;
      2 6;
      6 3;
      3 7;
      7 4 ];         % element connectivity matrix

nn=size(node,1);  % number of nodes
ndof=2*nn;        % number of dofs
ne=size(conn,1);  % number of elements

Ae=2.0*ones( 1, ne );       % element cross sectional areas
Ee=10e6*ones( 1, ne );      % element Young's modulus

fext=zeros(ndof,1); 

ifix=[1 2 8];               % the global dof ids of the fixed dofs
fext(12)=-10000             % the external force vector


% ----------------------- COMPUTATION SECTION -------------------------

fprintf('\n\n-----------------------------------------------\n');
fprintf(' 2D TRUSS FINITE ELEMENT CODE\n');

% assemble K
K=sparse(ndof,ndof);
for e=1:ne % loop over the elements
  
  conne=conn(e,:);  % element connectivity 
  
  % local stiffness matrix
  ke = ke_truss2d( node( conne,: ), Ae(e)*Ee(e) );  
  
  sctr=[ 2*conne(1)-1 2*conne(1) 2*conne(2)-1 2*conne(2) ];
  K(sctr,sctr) =  K(sctr,sctr) + ke;  % scatter ke into K
  
end

% solve the system
[d,R]=fesolve(K,fext,ifix);

% compute the strains and stresses and plot the deformed truss
clf
hold on
sfact=0.1*( max(max(node))-min(min(node)) )/max( [abs(d) ] );
for e=1:ne  
  
  n1=conn(e,1); % golbal node number for the first node in the element e
  n2=conn(e,2); % and for the second node in the element e
  
  X1= node(n1,:); % coordinate of the first node
  x1= X1+d(2*n1-1:2*n1)'; % deformed coordinate of the first node
  X2= node(n2,:); % coordinate of the second node
  x2= X2+d(2*n2-1:2*n2)'; % deformed coordinate of the second node
  
  l0e = norm( X2-X1 ); % undeformed length of the element
  le = norm( x2-x1 );  % deformed length of the element
  elong=le-l0e;
  
  strain(e) = elong/le;
  stress(e) = Ee(e)*strain(e);
  
  x1 = X1+sfact*d(2*n1-1:2*n1)'; % scaled deformed coord
  x2 = X2+sfact*d(2*n2-1:2*n2)'; 
  plot( [X1(1) X2(1)], [X1(2) X2(2)], 'r--' );
  plot( [x1(1) x2(1)], [x1(2) x2(2)], 'b-' );
  
end
axis equal

% print output
fprintf('\n  REACTION FORCES\n');
fprintf('%6d \t %6.3e\n',[ifix ;R']);

fprintf('\n  ELEMENT STRAIN AND STRESS\n');
fprintf('%6d \t %6.3e \t %6.3e\n',[1:ne ;strain; stress]);

fprintf('\n-----------------------------------------------\n');
