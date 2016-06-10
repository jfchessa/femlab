clear

% problem parameters
plateLength=100;
plateWidth=100;
plateThick=.125;

E=10e6;
nu=.3;
rho=.098/386;

% generate the mesh
nx=11; ny=11;
node=nodearray2d([0 0;plateLength,0;...
    plateLength,plateWidth;0,plateWidth],nx,ny);
element=genconn2d([1 2 nx+2 nx+1],nx-1,ny-1);
plot_mesh(node,element,'Quad4')
v=axis;
numElem=size(element,1);
numNode=size(node,1);
numDof=6*numNode;

% compute the stiffness and lumped mass matrices
Kmat=sparse(numDof,numDof);
Mmat=sparse(numDof,numDof);
for e=1:numElem
    econn = element(e,:);
    sctr = get_scatter( econn, 6 );
    ecoord = [node(econn,:),[0;0;0;0]];
    
    [ke, re, me] = kmat_mindlin_quad4(ecoord, plateThick, E, nu );
    Kmat(sctr,sctr) = Kmat(sctr,sctr) + ke;
    Mmat(sctr,sctr) = Mmat(sctr,sctr) + me*rho;
    
end

% solve the associated eigenvalue problem
[evect,lambda,iresult]=sptarn(Kmat,Mmat,1,1000,1);
nmodes=length(lambda);
for i=1:nmodes
    clf
    plot_mesh(node,element,'Quad4','k:')
    hold on
    sf=30;
    xcoord=node(:,1)+sf*evect(1:6:numDof,i);
    ycoord=node(:,2)+sf*evect(2:6:numDof,i);
    zcoord=sf*evect(3:6:numDof,i);
    natFreq=sqrt(lambda(i))/(2*pi);
    plot_mesh([node,zcoord],element,'Quad4','b-')
    title(['Mode ',num2str(i),', Natural Frequency=',num2str(natFreq),' [Hz]'])
    pause
end
    
