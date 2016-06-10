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

% Define the nodal coordinate matrix
node=[0 0;
    1 0;
    2 0;
    3 0;
    0 1;
    1 1;
    2 1;
    3 1 ];

% definte the element connectivity
conn=[ 1 2 6 5;
    2 3 7 6;
    3 4 8 7 ];

% define element properties
thk=0.125;
young=10e6;
poisson=0.33;

% define the fixed global dofs
ifix=[1 2 9 10];

% define the force vector
fext=zeros( size(node,1)*2, 1 );
fext(8)=-100;
fext(16)=-100;


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


