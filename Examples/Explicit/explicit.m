function status=explicit(jobname)
% -----------------------------------------------------------------------
%
%     Explicit finite element code
%
%
% -----------------------------------------------------------------------


% ---------------------  PREPROCESSOR  SECTION --------------------------
%clear
tic
status=1;

fprintf('\n----------------------------------------------------------------------------');
fprintf('\n|                                                                          |');
fprintf('\n|                                                                          |');
fprintf('\n|                                 FEMLAB CODE                              |');
fprintf('\n|                                                                          |');
fprintf('\n|      **************************************************************      |');
fprintf('\n|                                                                          |');
fprintf('\n|                               E X P L I C I T                            |');
fprintf('\n|                                                                          |');
fprintf('\n|                        F I N I T E      E L E M E N T                    |');
fprintf('\n|                                                                          |');
fprintf('\n|      **************************************************************      |');
fprintf('\n|                                                                          |');
fprintf('\n|       Written by Jack Chessa                                             |');
fprintf('\n|       jfchessa@utep.edu                                                  |');
fprintf('\n|       copyright 2010                                                     |');
fprintf('\n|                                                                          |');
fprintf('\n----------------------------------------------------------------------------');
fprintf('\n\n%8.2f: Reading input parameters\n',toc);

% define the mesh
nodeCoord=[0 0 0;
    5 0 0;
    10 0 0;
    0 1 0;
    5 1 0;
    10 1 0;
    0 0 1;
    5 0 1;
    10 0 1;
    0 1 1;
    5 1 1;
    10 1 1 ];

% define the elements
element(1) = struct('type','Hexa8',...
    'partid',1,...
    'conn',[1 2 5 4 7 8 11 10]);
element(2) = struct('type','Hexa8',...
    'partid',1,...
    'conn',[2 3 6 5 8 9 12 11]);

% define the parts
part(1) = struct('matid', 1, 'secid', 1);

% define the materials
material(1) = struct('type',@hypoelastic,'rho',.000245,'E',10e6,'nu',.33);

% and the section properties
section(1) = struct('type','solid');

% define the fixed global dofs
fixedNodes=[ 1   1   0.0 ;
    1   2   0.0 ;
    1   3   0.0 ;
    4   1   0.0 ;
    4   3   0.0 ;
    10   1   0.0 ;
    7   1   0.0 ];  % gnid ldof value

% timeinfo
tfin=.02;
cfl=0.5;

% define the nodal forces
iload=[2 1.25 0.0 0.0;    % global node and vector
    3 1.25 0.0 0.0 ];
sload=[1;1];             % scaling factor
tload=[0.0 0.0;
    .0001 1.0;
    tfin   1.0];       % temporal load curve

noutput=1;
nwrite=100;

% write fe geometry
%ensight_fegeometry('run.geom',nodeCoord,elemConn,'quad4');

% ------------------------  SOLVER SECTION  -----------------------------
nn=size(nodeCoord,1);
ndof=6*nn;
ne=size(element,1);

numHsv=1;
numQuadPts=8;

% compute fixed global dofs
ifix = 6*(fixedNodes(:,1)-1)+fixedNodes(:,2);
ival = fixedNodes(:,3);

% compute steady state fext
fext=zeros(ndof,1);
for i=1:size(iload,1)
    gnid=iload(i,1);
    fi=iload(i,2:4);
    fext( ((gnid-1)*6+1):((gnid-1)*6+3) ) = fi*sload(i);
end

% initialize system variables
mass=zeros(ndof,1);
fint=zeros(ndof,1);
dvect=zeros(ndof,1);
vvect=zeros(ndof,1);

% element data structures
hsve=zeros(numHsv,numQuadPts,ne);  % historical values
epspe=zeros(numQuadPts,ne);        % cumulative effective plastic strain
sige=zeros(6,numQuadPts,ne);       % effective plastic stress

fprintf('\n\n%8.2f: Problem summary',toc);
fprintf('\n\t\tNumber of nodes:%-6i',nn);
fprintf('\n\t\tNumber of elements:%-6i',ne);
fprintf('\n\t\tNumber of dofs:%-6i',ndof);
fprintf('\n\t\tFinal time:%-7.4e\n',tfin);
fprintf('\n\n%8.2f: Starting time integration\n',toc);

tn=0.0;
tstep=0;
dt=0;
dtlast=0;
wtimes=[];
while ( tn<=tfin )
    
    fint=zeros(ndof,1);
    tscale=1.0; %interp1(tload(:,1),tload(:,2),tn);
    
    for e=1:ne
        
        econn=element(e).conn;
        nne=length(econn);
        sctr=(econn-1)*6+1;
        sctr=reshape(ones(6,1)*sctr+(0:5)'*ones(1,nne),1,6*nne);
        ecoord=nodeCoord(econn,:);
        
        epartid = element(e).partid;
        ematid = part(epartid).matid;
        esecid = part(epartid).secid;
        etype=lower( element(e).type );
        
        switch lower( section(esecid).type )
            
            case 'solid'
                
               feval(['solid_',lower( element(e).type )]);
                
                
            case 'shell'
                
                
            case 'beam'
        
            otherwise
                disp(['Unknown element type ',section(esecid).type,...
                    ' in element ',num2str(e)]);
        end
        
%        [fe,me,dte,epspe(:,:,e), sige(:,:,e), hsve(:,:,e)]= ...
%            bt_shell( ecoord, vvect(sctr), elemThk, matLaw, ...
%            epspe(:,:,e), sige(:,:,e), hsve(:,:,e), rho0, matParam, ...
%            dtlast, numQuadPts );

        if ( tstep==0 )
            mass(sctr)=mass(sctr)+me;
        end
        
        fint(sctr)=fint(sctr)+fe;
        
        if ( e==1 || dte*cfl<dt )
            dt=dte*cfl;
        end
        
    end
    
    if ( mod(tstep,noutput)==0 )
        fprintf('\n%8.2f: Time step=%6i, sim. time=%7.4e, time step=%7.4e\n',toc,tstep,tn,dt);
        fprintf('          Element stress s11=%-8.2f s22=%-8.2f s12=%-8.2f\n',sige(1,1,e),sige(2,1,e),sige(4,1,e));
        if ( isnan(sige(1,1,e)) )
            break
        end
    end
    
    if ( mod(tstep,nwrite)==0 )
        fprintf('\n%8.2f: Time step=%6i, sim. time=%7.4e (writting data)\n',toc,tstep,tn);
        write_data(tstep/nwrite,sige,dvect);
        wtimes=[wtimes tn];
    end
    
    avect = (fext*tscale-fint)./mass;
    vvect = vvect + dt*avect;
    vvect(ifix)=0;
    dvect = dvect + dt*vvect;
    dvect(ifix)=0;
    
    tn=tn+dt;
    dtlast=dt;
    tstep=tstep+1;
    
end % temporal loop

fprintf('\n\n%8.2f: End of time integration',toc);
ensight_case( 'run','run.geom',wtimes,...
    {},{'disp'},{},...
    {'s11','s22','s33','s12','s23','s31'},{},{} );

fprintf('\n----------------------------------------------------------------------------');
fprintf('\n|                                                                          |');
fprintf('\n|   ****************      E N D    O F    R U N     ****************       |');
fprintf('\n|                                                                          |');
fprintf('\n----------------------------------------------------------------------------\n\n');

status =0;

end

function [ fe, me, dte, epspe, sige, hsve ] = solid_hexa8( ecoord, ve, matlaw, eparam )

    fe=0;
    me=0;
    dte=0;
    epse=0;
    sige=0;
    hsve=0;

end

