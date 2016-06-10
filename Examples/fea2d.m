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
%----------------------------------------------------------------------


fprintf('\n\n-----------------------------------------------\n');
fprintf('|    2D PLANE STRESS FINITE ELEMENT CODE      |\n');
fprintf('-----------------------------------------------\n');


nn=size(node,1);  % number of nodes
ndof=2*nn;        % number of dofs
ne=length(element);  % number of elements

% ----------------------- COMPUTATION SECTION -------------------------

% assemble K
fprintf(' ASSEMBLING STIFFNESS MATRIX\n');

c1=young/(1-poisson^2);  % plne stress material stiffness matrix
c2=poisson*c1;
c3=0.5*(1-poisson)*c1;
C=[c1 c2 0;c2 c1 0; 0 0 c3];

for e=1:ne

  pid=element(e).pid;
  mid=part(pid).mid;
  sid=part(pid).sid;
  
  young=material(mid).young_mod;
  nu=material(mid).poisson_ratio;
  
  conne=element(e).conn;
  nne=length(conne);
  
  coord=node(conne,1:2);

  sctr(1:2:2*nne)=dof_map(conne,1);
  sctr(2:2:2*nne)=dof_map(conne,2);
  
  switch element(e).type

    % plane stress elements
    case 'Tria3'
      thk=section(sid).thk;
      ke=kmat_tria3( coord, young, nu, thk );

    case 'Tria6'
      thk=section(sid).thk;
      ke=kmat_tria6( coord, young, nu, thk );

    case 'Quad4'
      thk=section(sid).thk;
      ke=kmat_quad4( coord, young, nu, thk );

    case 'Quad8'
      thk=section(sid).thk;
      ke=kmat_quad8( coord, young, nu, thk );

    % structural elements
    case 'Line2'
      if ( section(sid).type=='truss' )
        A=section(sid).area;
        ke=kmat_truss2d( coord, A*young );

      else ( section(sid).type=='beam' )
        A=section(sid).area;
        I=section(sid).ixx;
        ke=kmat_beam2d( coord, I*young, A*young );
        sctr(1:3:3*nne)=dof_map(conne,1);
        sctr(2:3:3*nne)=dof_map(conne,2);
        sctr(3:3:3*nne)=dof_map(conne,3);
      end

    otherwise
      disp(['Error element',num2str(e),....
        ' unsupported element type ',element(e).type]);

  end

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
