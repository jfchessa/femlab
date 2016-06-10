%
%  explicit1d.m
%
%
%
clear

barLength=10;
rho=1.0;
young=500.0;
area=1.0;

tfin=1.0;
CFL=0.9;

stress0=10;

numElem=500;

% -----------------
numNode=numElem+1;
node=linspace(0,barLength,numNode);
conn=[1:numNode-1;2:numNode]';
numDof=numNode;

% Form the stiffness matrix
KMat=sparse(numDof,numDof);
for e = 1:numElem
    econn=conn(e,:);
    le=abs( node(econn(2))-node(econn(1)) );
    ke=area*young/le*[1 -1;-1 1];
    sctr=econn;
    KMat(sctr,sctr) = KMat(sctr,sctr) + ke;
end

% and the mass matrix
MMat=sparse(numDof,numDof);
for e = 1:numElem
    econn=conn(e,:);
    le=abs( node(econn(2))-node(econn(1)) );
    me=0.5*rho*area*le*[1 0;0 1];
    sctr=econn;
    MMat(sctr,sctr) = MMat(sctr,sctr) + me;
end

% and the stable time step
dt=10e10;
waveSpeed=sqrt(young/rho);
for e = 1:numElem
    econn=conn(e,:);
    le=abs( node(econn(2))-node(econn(1)) );
    dte=CFL*le/waveSpeed;
    if dte<dt
        dt=dte;
    end
end


% initial conditions
dvect=zeros(numDof,1);
vvect=zeros(numDof,1);
fext=zeros(numDof,1);
fint=zeros(numDof,1);
avect=zeros(numDof,1);
tn=0.0;

clf
while tn<tfin
    
    fext(1)=stress0*area;
    
    fint=KMat*dvect;
    avect = MMat\(fext-fint);
    
    vvect = vvect + dt*avect;
    dvect = dvect + dt*vvect;
    tn = tn + dt;
    
    plot(node,vvect)
    grid on
    vmax=stress0/sqrt(young*rho);
    axis([0,barLength,-vmax,3*vmax])
    xlabel('Distance along bar')
    ylabel('Velocity')
    pause(0.1)
    
end