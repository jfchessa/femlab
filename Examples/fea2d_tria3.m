%----------------------------------------------------------------------
%
% fea2d.m
%
% A two dimensional finite element code 
%
% Written by Jack Chessa jfchessa@utep.edu
% for CE 5307 Theory of Finite Element Analysis
% November 6, 2007
%
%---------------------------------------------------------------------
% 
%  Input Variables (must be defined)
%
%       node - nodal coordinate matrix (2xnn matrix)
%       conn - element connectivity matrix (2xnumelem matrix)
%       area - element cross sectional area (column vector)
%       young - element Young's modulus (column vector)
%
%       ifix - fixed global dofs
%       iforce - global dofs where point forces are applied (1st col)
%                and value of nodal forces (2nd col)
%
%       defScale - optional parameters that scales teh displacement on 
%                  output plot
%
%---------------------------------------------------------------------
clear
%           S T A R T    O F   I N P U T     S E C T I O N
node=[0.0 0.0;
      0.5 0.0;
      1.0 0.0;
      0.0 0.5;
      0.5 0.5;
      1.0 0.5;
      0.0 1.0;
      0.45 1.05;
      0.9 1.1;
      0.0 1.5;
      0.5 1.5;
      1.0 1.5;
     1.15 0.85;
      1.5 1.0;
      2.0 1.0;
      1.5 1.5;
      2.0 1.5 ];
  
element=[1 2 5; 
    2 3 6;
    4 5 8;
    5 6 9;
    7 8 11;
    8 9 12;
    1 5 4;
    2 6 5;
    4 8 7;
    5 9 8;
    7 11 10;
    8 12 11;
    6 13 9;
    13 14 9;
    9 14 16;
    14 15 17;
    9 16 12;
    14 17 16 ];  

young=10e6;
poisson=0.33;

thk=1.0;

iforce = [ 15*2 500; 
           17*2 500];
       
ifix = [1 2  2*2 2*3 2*4-1 2*7-1 2*10-1];

defScale = 100;
%
%             E N D     O F   I N P U T     S E C T I O N
%---------------------------------------------------------------------
%
% You should not need to touch anything below
%
fprintf('\n\n-----------------------------------------------\n');
fprintf('|    2D PLANE STRESS FINITE ELEMENT CODE      |\n');
fprintf('-----------------------------------------------\n');

nn=size(node,1);  % number of nodes
ndof=2*nn;        % number of dofs
ne=size(element,1);  % number of elements

% ----------------------- COMPUTATION SECTION -------------------------

% assemble K
fprintf(' ASSEMBLING STIFFNESS MATRIX\n');

c1=young/(1-poisson^2);  % plne stress material stiffness matrix
c2=poisson*c1;
c3=0.5*(1-poisson)*c1;
C=cmat_mat1(young,poisson,'PSTRESS');

K=sparse(ndof,ndof);
for e=1:ne

  conne=element(e,:);
 
  sctr(1:2:6)=2*conne-1;
  sctr(2:2:6)=2*conne;

  [B,A]=bmat_tria3( node(conne,:) );
  ke=B'*C*B*A*thk;
  
  K(sctr,sctr) = K(sctr,sctr) + ke;	

end

% compute the external force
fext=zeros(ndof,1);
fext(iforce(:,1))=iforce(:,2);

% solve the system
fprintf(' ASSEMBLING STIFFNESS MATRIX\n');
[d,freac]=fesolve(K,fext,ifix);

% compute the strains and stresses and plot results
fprintf(' CALCULATE THE STRESSES\n');
stress=zeros(4,ne);  % matrix of element stresses and strains
strain=zeros(3,ne);  %

for e=1:ne

  conne=element(e,:);
  B=bmat_tria3( node(conne,:) );  
  sctr(1:2:6)=2*conne-1;
  sctr(2:2:6)=2*conne;
  strain(:,e)=B*d(sctr);
  stress(1:3,e)=C*strain(:,e); 
  
  ps=principal_val(stress(1:3,e));  % principal streses
  stress(4,e)=sqrt( 0.5*( (ps(2)-ps(1))^2 + (ps(3)-ps(2))^2 + ...
    (ps(1)-ps(3))^2 ) );   % Mises stress
  
  
end

% ---------------------- POST PROCESSING SECTION -------------------------
% plot the stresses
if ( ~exist('defScale') )
    defScale = 10;
end
x = node + defScale*[d(1:2:ndof) d(2:2:ndof)];
fprintf(' PLOTTING RESULTS\n');
clf
hold on
trisurf(element,x(:,1),x(:,2),zeros(nn,1),stress(4,:),'EdgeColor','cyan')
trimesh(element,node(:,1),node(:,2),zeros(nn,1),'FaceColor','none','EdgeColor','black','LineStyle','--')
title('PLOT OF MISES STRESS'); 
colorbar
view(2)
axis equal

% write the results to an ensight file (you can  read this with
% paraview,  www.paraview.org)
fprintf(' WRITING RESULTS\n');
ensight_fegeometry('fea2d.geom',node,element,'Tria3');
ensight_field('fea2d0000.s11',nodal_avg(stress(1,:),element,node));
ensight_field('fea2d0000.s22',nodal_avg(stress(2,:),element,node));
ensight_field('fea2d0000.s12',nodal_avg(stress(3,:),element,node));
ensight_field('fea2d0000.svm',nodal_avg(stress(4,:),element,node));
ensight_field('fea2d0000.e11',nodal_avg(strain(1,:),element,node));
ensight_field('fea2d0000.e22',nodal_avg(strain(2,:),element,node));
ensight_field('fea2d0000.e12',nodal_avg(strain(3,:),element,node));
ensight_case('fea2d','fea2d.geom',0,...
  {'s11','s22','s12','svm','e11','e22','e12'});

fprintf('\n-----------------------------------------------\n');
fprintf(  '|                END OF PROGRAM               |');
fprintf('\n-----------------------------------------------\n');
