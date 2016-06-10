%----------------------------------------------------------------------
%
% sol101.m
%
% A linear elastic finite element code 
%
%
%
%----------------------------------------------------------------------

%
% This code mimics some of the capabilities of a NASTRAN SOL101
%

fprintf('\n\n-----------------------------------------------\n');
fprintf('|                FEMLAB SOL 101               |\n');
fprintf('-----------------------------------------------\n');

% Define the basic data structures
%
%  nodeCoord = The nodal coordinate matrix
%
%  rodElem = An array of element structures
%  barElem = An array of element structures
%  shellElem = An array of element structures
%  solidElem = An array of element structures
%
%  matData =
% 
%  rodProp = 
%  barProp = 
%  shellProp = 
%  solidProp = 
%
%  forceData = 
%  spcData = 
%
%  loadCase = 


% read the NASTRAN input section
%  
%  All nodes are renumbered to be zero offset and consecutive
%
%
fprintf(' READING INPUT DATA\n');

% some fake in data
nodeCoord=[ 0.0, 0.0, 0.0;
            1.0, 0.0, 0.0;
            1.0, 1.0, 0.0;
            0.0, 1.0, 0.0;
            0.0, 0.0, 2.0;
            1.0, 0.0, 2.0;
            1.0, 1.0, 2.0;
            0.0, 1.0, 2.0 ];
        
shellElem(1).type='Quad4';
shellElem(1).pid=1;
shellElem(1).conn=[1 2 6 5];

shellElem(2).type='Quad4';
shellElem(2).pid=1;
shellElem(2).conn=[3 4 8 7];

barElem(1).type='Line2';
barElem(1).pid=1;
barElem(1).conn=[2 3];
barElem(1).v=[1 0 0];

rodElem(1).type='Line2';
rodElem(1).pid=1;
rodElem(1).conn=[3 7];

matData(1).type=1;
matData(1).young=10e6;
matData(1).poisson=0.3;

matData(2).type=1;
matData(2).young=30e6;
matData(2).poisson=0.28;

rodProp(1).mid=2;
rodProp(1).area=1.0;
rodProp(1).J=.25;
rodProp(1).c=.125;

barProp(1).mid=2;
barProp(1).area=1.0;
barProp(1).I11=.001;
barProp(1).I22=.05;
barProp(1).J=.051;
rodProp(1).c1=[-.06 0];
rodProp(1).c2=[ .06 0];
rodProp(1).c3=[ 0 .06];
rodProp(1).c4=[ 0 -.06];

shellProp(1).mid=1;
shellProp(1).thk=.125;

forceData(1).nid=[7 8]
forceData(1).value=[0 100 0 0 0 0; 0 50 0 0 0 0];

forceData(2).nid=[8]
forceData(2).value=[0 -100 0 0 0 0];

spcData(1).nid=[1 2 3 4];
spcData(1).dofs=['1','2','3','123'];

loadCase(1).spcid=1;
loadCase(1).loadid=1;

loadCase(2).spcid=1;
loadCase(2).loadid=2;


% ----------------------- COMPUTATION SECTION -------------------------

% setup dof map
[dofMap,ngdof]=set_dofmap(solidlElem,shellElem,barElem,rodElem);

% compute the load vector
numLC=length(loadCase);
fext=zeros(ngdof,numLC);
fext = compute_force( loadCase, dofMap, forceData, fext );

% compute the fixed global dofs
[ifix,ival]=get_spcs( loadCase, dofMap, spcData );

% assemble K
fprintf(' ASSEMBLING STIFFNESS MATRIX\n');
K=sparse(ngdof,ngdof);
K=compute_solid( nodeCoord, dofMap, solidlElem, solidProp, matData, K );
K=compute_shell( nodeCoord, dofMap, shellElem, shellProp, matData, K );
K=compute_bar( nodeCoord, dofMap, barElem, barProp, matData, K );
K=compute_rod( nodeCoord, dofMap, rodElem, rodProp, matData, K );

% solve the system
fprintf(' SOLVEING THE SYSTEM EQUATIONS\n');
[d,freac]=fesolve(K,fext,ifix,ival);

% compute the strains and stresses and plot results
fprintf(' CALCULATE THE STRESSES\n');


% ---------------------- POST PROCESSING SECTION -------------------------

% write the results to an ensight file (you can  read this with
% paraview,  www.paraview.org)
fprintf(' WRITING RESULTS\n');


fprintf('\n-----------------------------------------------\n');
fprintf(  '|                END OF PROGRAM               |');
fprintf('\n-----------------------------------------------\n');
