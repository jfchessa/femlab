% femQ4.m
%
% A simple static 2d finite elment program that uses Q4 elements
%
clear

%----------------------------------------------------------------------
%
% INPUT SECTION
%
%

L=10;
H=4;

he=0.5;

% Define the nodal coordinate matrix
nex=ceil(L/he);
ney=ceil(H/he);
node=nodearray2d( [0 0;L 0;L H; 0 H], nex+1, ney+1 );
nn=size(node,1);

% definte the element connectivity
conn=genconn2d([1 2 nex+3 nex+2],nex,ney,1,2);
plot_mesh(node,conn,'Quad4')

% define element properties
thk=0.125;
young=10e6;
poisson=0.33;

% define the fixed global dofs
lhsNodes=1:(nex+1):nn;
ifix=[2*lhsNodes-1,2];

% define the force vector
rhsNodes=(nex+1):(nex+1):nn;
fext=zeros( 2*nn, 1 );
P=-100;
fext(2*rhsNodes)=P/(H/ney);
fext(2*rhsNodes(1))=0.5*P/(H/ney);
fext(2*nn)=0.5*P/(H/ney);


%----------------------------------------------------------------------
%
% COMPUTATIONAL SECTION
%
%

numNode=size(node,1);
numElem=size(conn,1);
numDof=2*numNode;

% compute the element stiffness matrix
K=sparse(numDof,numDof);
for e=1:numElem
    
    coord=node(conn(e,:),:);  % matrix of the element nodal coordinates
    
    sctr(1:2:7) = 2*conn(e,:)-1;  % scatter/gather vector
    sctr(2:2:8) = 2*conn(e,:);
    
    ke = kmat_quad4(coord, young, poisson, thk);  % conpute the element 
                                                  % stiffness matrix
    
    K(sctr,sctr)=K(sctr,sctr)+ke;   % scatter the element ke into K
                                                  
end

% solve the system
[d,freac]=fesolve(K,fext,ifix);

%----------------------------------------------------------------------
%
% POST-PROCESSING SECTION
%
%

% write the nodal displacements
fprintf('\n\n  nid   x-displacement  y-displacement\n')
fprintf('--------------------------------------\n')
for i=1:numNode
    fprintf('%6i    %12.6g    %12.6g \n',i,d(2*i-1),d(2*i));
end


% compute and write the element stresses
estress=zeros(numElem,4);
fprintf('\n\n   eid          S11          S22          S12          Svm\n')
fprintf('----------------------------------------------------------\n')
for e=1:numElem
    
    coord=node(conn(e,:),:);  % matrix of the element nodal coordinates
    
    sctr(1:2:7) = 2*conn(e,:)-1;  % scatter/gather vector
    sctr(2:2:8) = 2*conn(e,:);
    
    B = bmat_quad4(coord,[0 0]);  % compute the element B matrix
    strain=B*d(sctr);
    
    C = cmat_mat1( young, poisson, 'pstress' );
    stress=C*strain;
    estress(e,1:3)=stress';
    estress(e,4)=mises_val(stress);
    
    fprintf('%6i    %9.3g    %9.3g    %9.3g    %9.3g     \n',...
        e, stress(1), stress(2), stress(3), estress(e,4) );
    
end



% plot the results
svm=nodal_avg(estress(:,4),conn,node);  
plot_fefield(node,conn,svm,'Quad4')
title('Plot of Mises Stress')
colorbar


