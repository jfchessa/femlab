%----------------------------------------------------------------------
%
% heat_conduction.m
%
% A two dimensional heat conduction finite element code 
%
% Written by Jack Chessa jfchessa@utep.edu
% for CE 5307 Theory of Finite Element Analysis
% November 6, 2007
%
%----------------------------------------------------------------------


fprintf('\n\n-----------------------------------------------\n');
fprintf('|    2D HEAT CONDUCTION FINITE ELEMENT CODE      |\n');
fprintf('-----------------------------------------------\n');


nn=size(node,1);     % number of nodes
ndof=nn;             % number of dofs
ne=size(element,1);  % number of elements

% ----------------------- COMPUTATION SECTION -------------------------

% assemble K
fprintf(' ASSEMBLING STIFFNESS MATRIX\n');

Conductivity=100;
HeatCapacitance=1;
Density=1;
alpha=Conductivity/HeatCapacitance/Density;

K=sparse(ndof,ndof);
for e=1:ne

  conne=element(e,:);
  sctr=conne;

  [B,A]=bmat_tria3( node(conne,:) );
  ke=B'*alpha*B*A*thk;
  
  K(sctr,sctr) = K(sctr,sctr) + ke;	

end

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
fprintf(' PLOTTING RESULTS\n');
clf
trisurf(element,node(:,1),node(:,2),zeros(nn,1),stress(4,:))
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
