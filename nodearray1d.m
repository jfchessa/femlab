function a=nodearray1d(p1,p2,n,stype,sparam)

% NODES=NODEARRAY1D(P1,P2,STYPE,SPARAM) 
%
% Creates a mesh along a line
%
%  P1, P2 the points between which the nodes are to be generated.
%  N - the number of nodes between the points
%  STYPE - scaling type 'LINEAR', 'POWER', BELLCURVE'
%  SPARAM - parameter for the scaling

xi=linspace(0,1,n)';

if ( nargin<4 )
    stype='NONE';
end
if ( nargin<5 )
    sparam=1;
end

if ( strcmp(lower(stype),'linear') )
    rv=sparam^(1/(n-2));   % get scaling
    xi(1)=0;
    d=1;
    for i=2:n
        xi(i)=xi(i-1)+d;
        d=d/rv;
    end
    xi=xi/xi(n);
elseif ( strcmp(lower(stype),'power') )
    xi=xi.^sparam;
elseif ( strcmp(lower(stype),'bellcurve') )
    bc=0.5*(tanh(sparam*(xi-0.5))+1);
    bc=bc-bc(1);
    bc=bc/bc(n);
    xi=bc;
end

a=(ones(n,1)-xi)*p1+xi*p2;  

