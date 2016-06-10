%----------------------------------------------------------------------
%
% fea1d.m
%
% A finite element code for one dimensional bar analysis
%
% Written by Jack Chessa jfchessa@utep.edu
% for CE 5307 Theory of Finite Element Analysis
% September 6, 2007
%
%----------------------------------------------------------------------
clear

% ------------------------ DATA INPUT SECTION -------------------------

node=[0.0; 30; 75; 115];  % node coordinate matrix
conn=[1 2;2 3;3 4];         % element connectivity matrix

Ae=[3.0 3.0 3.0];           % element cross sectional areas
Ee=30e6*[1 1 1];            % element Young's modulus

ifix=[1];                     % the global node ids of the fixed nodes
fext=[ 0.0; 81000; -66000; 30000 ]; % the external force vector

% ----------------------- COMPUTATION SECTION -------------------------

fprintf('\n\n-----------------------------------------------\n');
fprintf(' 1D FINITE ELEMENT CODE\n');
nn=size(node,1);
ne=size(conn,1);

% assemble K
K=sparse(nn,nn);
for e=1:ne % loop over the elements
  
  n1=conn(e,1); % golbal node number for the first node in the element e
  n2=conn(e,2); % and for the second node in the element e
  
  x1= node(n1); % coordinate of the first node
  x2= node(n2); % coordinate of the second node
  
  le = abs( x2-x1 ); % length of the element
  
  ke = Ae(e)*Ee(e)/le*[1.0 -1.0;-1.0 1.0];  % local stiffness matrix
  
  sctr=conn(e,:);  % local to global scatter vector
  K(sctr,sctr) =  K(sctr,sctr) + ke;  % scatter ke into K
  
end

% solve the system
[d,R]=fesolve(K,fext,ifix);

% compute the strains and stresses
for e=1:ne  
  
  n1=conn(e,1); % golbal node number for the first node in the element e
  n2=conn(e,2); % and for the second node in the element e
  
  x1= node(n1); % coordinate of the first node
  x2= node(n2); % coordinate of the second node
  
  le = abs( x2-x1 ); % length of the element
  elong = d(n2)-d(n1); % elongation
  strain(e) = elong/le;
  stress(e) = Ee(e)*strain(e);
end

% print output
fprintf('\n  NODAL DISPLACEMENTS\n');
fprintf('%6d \t %6.3e\n',[1:nn ;d']);

fprintf('\n  REACTION FORCES\n');
fprintf('%6d \t %6.3e\n',[ifix ;R']);

fprintf('\n  ELEMENT STRAIN AND STRESS\n');
fprintf('%6d \t %6.3e \t %6.3e\n',[1:ne ;strain; stress]);

fprintf('\n-----------------------------------------------\n');
