%----------------------------------------------------------------------
%
% fea2d_axisym.m
%
% A two dimensional axisymmetric finite element code 
%
% Written by Jack Chessa jfchessa@utep.edu
%
%----------------------------------------------------------------------
clear

%----------------------------------------------------------------------
%
%   INPUT SECTION
%
Rout=15;
Rin=0.0;
HeightSubstrait=2;
ThicknessOxide=1;
WaveAmp=0.0;
WaveLength=1;

ElementSize=0.125;

youngSubstrait=10e6;
poissonSubstrait=0.33;
alphaSubstrait=.0002;

youngOxide=10e6;
poissonOxide=0.33;
alphaOxide=.0001;

dT=200;

%----------------------------------------------------------------------
%
% Should not need to touch anything after this
%
Length=Rout-Rin;
nr=ceil(Length/ElementSize)+1;
nz1=ceil(HeightSubstrait/ElementSize)+1;
nz2=ceil(ThicknessOxide/ElementSize)+1;
hzOxide=ThicknessOxide/(nz2-1);


node=nodearray2d( [ Rin  0;
                    Rout 0; 
                    Rout HeightSubstrait; 
                    Rin  HeightSubstrait], nr, nz1 );
                
node = [ node;   
    nodearray2d( [ Rin  HeightSubstrait+hzOxide;
                   Rout HeightSubstrait+hzOxide; 
                   Rout HeightSubstrait+ThicknessOxide; 
                   Rin  HeightSubstrait+ThicknessOxide], nr, nz2-1 )   ];
          
                
element=genconn2d([1 2 nr+2 nr+1],nr-1,nz1+nz2-2,1,2);
OxideElemStart=(nr-1)*(nz1-1)+1;

nn=size(node,1);  
ndof=2*nn;         
ne=size(element,1);   

issub=find(node(:,2)<HeightSubstrait);
node(issub,2) = node(issub,2).*(1 + ...
    WaveAmp/HeightSubstrait*sin(2*pi*node(issub,1)/WaveLength) );
isoxide=max(issub)+1:nn;
node(isoxide,2) = node(isoxide,2) + ...
    WaveAmp*sin(2*pi*node(isoxide,1)/WaveLength);

nid_left=1:nr:nn;
nid_right=nr:nr:nn;

%ifix=[ 2:2:2*nr 1:2*nr:ndof ];
%ifix=[ 2:2:2*nr ];
ifix=2;
fext=zeros(ndof,1);

xdofs=2*nid_right-1;
xcorners=xdofs([1 length(xdofs)]);
fext(xdofs)=-100*hzOxide*(2*pi*Rout);
fext(xcorners)=fext(xcorners)/2;

%xdofs=2*nid_left-1;
%xcorners=xdofs([1 length(xdofs)]);
%fext(xdofs)=100*hzOxide*(2*pi*Rin);
%fext(xcorners)=fext(xcorners)/2;

thk=1.0;

clf
%plot_mesh(node,element,'quad4');
%break

fprintf('\n\n-----------------------------------------------\n');
fprintf(    '|    2D AXISYMMETRIC FINITE ELEMENT CODE      |\n');
fprintf(    '-----------------------------------------------\n');

% ----------------------- COMPUTATION SECTION -------------------------

% assemble K
fprintf(' ASSEMBLING STIFFNESS MATRIX\n');


K=sparse(ndof,ndof);
[ qpts, qwts ] = element_quadrature( 'quad4', 2 );
%[ zpts, zwts ] = element_quadrature( 'line2', 2 );
%[ rpts, rwts ] = element_quadrature( 'line2', 12 );
%[ qpts, qwts ] = quadrature_compound( rpts, rwts, zpts, zwts );
% multiply weights by 2 pi
qwts=2*pi*qwts;

C=cmat_mat1(youngSubstrait,poissonSubstrait,'AXISYMMETRIC');
alpha=alphaSubstrait;
for e=1:ne

  conne=element(e,:);
 
  sctr(1:2:8)=2*conne-1;
  sctr(2:2:8)=2*conne;
  
  if ( e==OxideElemStart )
    C=cmat_mat1(youngOxide,poissonOxide,'AXISYMMETRIC');
    alpha=alphaOxide;
  end
  
  for q=1:length(qwts)
      
    [B,rpt,jac,N]=bmat_axisymmetric( node(conne,:), qpts(q,:), 'quad4' );
    Tpt=dT;
    thermalStrain=alpha*[1;1;1;0]*Tpt;
    
    K(sctr,sctr) = K(sctr,sctr) + B'*C*B*jac*qwts(q)*rpt;
    fext(sctr) = fext(sctr) + B'*C*thermalStrain*jac*qwts(q)*rpt;
  
  end

end

% solve the system
fprintf(' SOLVING THE SYSTEM EQUATIONS\n');
[d,freac]=fesolve(K,fext,ifix);

% compute the strains and stresses and plot results
fprintf(' CALCULATE THE STRESSES\n');
stress=zeros(5,ne);  % matrix of element stresses and strains
strain=zeros(4,ne);  %

C=cmat_mat1(youngSubstrait,poissonSubstrait,'AXISYMMETRIC');
alpha=alphaSubstrait;
for e=1:ne
    
  conne=element(e,:);
  sctr(1:2:8)=2*conne-1;
  sctr(2:2:8)=2*conne;
  
  if ( e==OxideElemStart )
    C=cmat_mat1(youngOxide,poissonOxide,'AXISYMMETRIC');
    alpha=alphaOxide;
  end
  
  Tpt=dT;
  thermalStrain=alpha*[1;1;1;0]*Tpt;
      
  xi = [0,0];
  [B,rpt,A]=bmat_axisymmetric( node(conne,:), xi, 'quad4' );
  
  strain(:,e)=B*d(sctr)-thermalStrain;
  stress(1:4,e)=C*strain(:,e); 
  
  sv=[ stress(1:3,e); 0; stress(4,e); 0]; 
  ps=principal_val(sv);  % principal streses
  stress(5,e)=sqrt( 0.5*( (ps(2)-ps(1))^2 + (ps(3)-ps(2))^2 + ...
    (ps(1)-ps(3))^2 ) );   % Mises stress
  
  
end

% ---------------------- POST PROCESSING SECTION -------------------------
% plot the stresses
fprintf(' PLOTTING RESULTS\n');
clf
trisurf(element,node(:,1)+d(1:2:ndof),node(:,2)+d(2:2:ndof),zeros(nn,1),stress(5,:))
title('PLOT OF MISES STRESS'); 
colorbar
view(2)
axis equal

% write the results to an ensight file (you can  read this with
% paraview,  www.paraview.org)
fprintf(' WRITING RESULTS\n');
ensight_fegeometry('fea2d.geom',node,element,'Quad4');
ensight_field('fea2d0000.dsp', [d(1:2:ndof) d(2:2:ndof)] );
ensight_field('fea2d0000.srr',nodal_avg(stress(1,:),element,node));
ensight_field('fea2d0000.stt',nodal_avg(stress(2,:),element,node));
ensight_field('fea2d0000.szz',nodal_avg(stress(3,:),element,node));
ensight_field('fea2d0000.srz',nodal_avg(stress(4,:),element,node));
ensight_field('fea2d0000.svm',nodal_avg(stress(5,:),element,node));
ensight_field('fea2d0000.err',nodal_avg(strain(1,:),element,node));
ensight_field('fea2d0000.ett',nodal_avg(strain(2,:),element,node));
ensight_field('fea2d0000.ezz',nodal_avg(strain(3,:),element,node));
ensight_field('fea2d0000.erz',nodal_avg(strain(4,:),element,node));
ensight_case('fea2d','fea2d.geom',0,...
  {'srr','stt','szz','srz','svm','err','ett','ezz','erz'},...
  {'dsp'});

fprintf('\n-----------------------------------------------\n');
fprintf(  '|                END OF PROGRAM               |');
fprintf('\n-----------------------------------------------\n');
