function [ke, re, me] = kmat_mindlin_quad4(gcoord, h, e, nu, qvect,...
    dir, rotzStiffFactor)
% [ke, re] = kmat_mindlin_quad4(gcoord, h, e, nu, qvect, dir, rotzStiffFactor)
%
% Function to compute element stiffness and equivalent load vector for a
% four node quadralateral Mindlin quad4 shell element
%
% 	gcoord = element global coordinates
% 	h = thickness
% 	e = Young's modulus, nu = Poisson's ratio
% 	qvect = distributed load along three coordinate directions.
% 	dir = load direction. If dir = "Local" the distributed loads are
% 				specified in the element local coordinates, otherwise they are assumed
% 				to be in the global coordinates.  'Global' is the default
% 	rotzStiffFactor = normal rotational stiffness factor
%
% 	ke =  element stiffness, re =  equivalent load vector.
%   me =  lumped mass matrix/density
%
%
% [ke, re, me] = kmat_mindlin_quad4( gcoord, h, e, nu )
% [ke, re, me] = kmat_mindlin_quad4( gcoord, h, e, nu, qvect )
% [ke, re, me] = kmat_mindlin_quad4( gcoord, h, e, nu, qvect, dir )
%
%
%

if ( nargin < 5 )
	q1=0; q2=0; q3=0;
else
	q1=qvect(1); q2=qvect(2); q3=qvect(3); 
end

if ( nargin< 6 )
	dir = 'Global';
end

if ( nargin<7 )
	rotzStiffFactor=1;
end

elemArea = 0;
H = ShellDirectionCosines(gcoord);  % direction cosines of the element
coord=gcoord;

for i=1:size(gcoord,1)
    coord(i,:) = (H*gcoord(i,:)')';
end

qlocal = [q1; q2; q3];
if strcmp(dir,'Global')
	qlocal = H*qlocal;
end

qs = qlocal(1); qt = qlocal(2); qr=qlocal(3);
T = zeros(24);
for i=1:3:24
    T(i:i+2,i:i+2)=H;
end

% elastic constants
g = e/(2*(1 + nu)); d = e*h^3/(12*(1 - nu^2));
dp = e*h /(1 - nu^2);
c = [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu)/2];

% Gauss point locations and weights
pt = 1/sqrt(3);
gBendLocs = [-pt, -pt; pt, -pt; -pt, pt; pt, pt];
gBendWts = [1, 1, 1, 1];
pt = 0;
gShearLocs = [[pt, pt]];
gShearWts = [4];
dof = 24;
ke=zeros(dof); re=zeros(dof,1);
k1=zeros(8); k2=zeros(12);
for i=1:size(gBendLocs,1)
    s = gBendLocs(i, 1); t = gBendLocs(i, 2); w = gBendWts(i);
    n = [(1/4)*(1 - s)*(1 - t), (1/4)*(s + 1)*(1 -t), ...
        (1/4)*(s + 1)*(t + 1), (1/4)*(1 - s)*(t + 1)];
    dns=[(-1 + t)/4, (1 - t)/4, (1 + t)/4, (-1 - t)/4];
    dnt=[(-1 + s)/4, (-1 - s)/4, (1 + s)/4, (1 - s)/4];
    x = n*coord(:,1); y = n*coord(:,2);
    dxs = dns*coord(:,1); dxt = dnt*coord(:,1);
    dys = dns*coord(:,2); dyt = dnt*coord(:,2);
    J = [dxs, dxt; dys, dyt]; detJ = det(J);
    bx = (J(2, 2)*dns - J(2, 1)*dnt)/detJ;
    by = (-J(1, 2)*dns + J(1, 1)*dnt)/detJ;
    b = [bx(1),0,bx(2),0,bx(3),0,bx(4),0;
        0,by(1),0,by(2),0,by(3),0,by(4);
        by(1),bx(1),by(2),bx(2),by(3),bx(3),by(4),bx(4)];
    k1 = k1 + dp*detJ*w*b'*c*b;
    b = [0,0,-bx(1),0,0,-bx(2),0,0,-bx(3),0,0,-bx(4);
        0,by(1),0,0,by(2),0,0,by(3),0,0,by(4),0;
        0,bx(1),-by(1),0,bx(2),-by(2),0,bx(3),-by(3),0,bx(4),-by(4)];
    k2 = k2 + d*detJ*w*b'*c*b;
    re = re + detJ*w*[qs*n(1),qt*n(1),qr*n(1),0,0,0,...
        qs*n(2),qt*n(2),qr*n(2),0,0,0,...
        qs*n(3),qt*n(3),qr*n(3),0,0,0,...
        qs*n(4),qt*n(4),qr*n(4),0,0,0]';
elemArea = elemArea + detJ*w;
end

for i=1:size(gShearLocs,1)
    s = gShearLocs(i, 1); t = gShearLocs(i, 2); w = gShearWts(i);
    n = [(1/4)*(1 - s)*(1 - t), (1/4)*(s + 1)*(1 -t), ...
        (1/4)*(s + 1)*(t + 1), (1/4)*(1 - s)*(t + 1)];
    dns=[(-1 + t)/4, (1 - t)/4, (1 + t)/4, (-1 - t)/4];
    dnt=[(-1 + s)/4, (-1 - s)/4, (1 + s)/4, (1 - s)/4];
    x = n*coord(:,1); y = n*coord(:,2);
    dxs = dns*coord(:,1); dxt = dnt*coord(:,1);
    dys = dns*coord(:,2); dyt = dnt*coord(:,2);
    J = [dxs, dxt; dys, dyt]; detJ = det(J);
    bx = (J(2, 2)*dns - J(2, 1)*dnt)/detJ;
    by = (-J(1, 2)*dns + J(1, 1)*dnt)/detJ;
    b = [bx(1),0,n(1),bx(2),0,n(2),bx(3),0,n(3),bx(4),0,n(4);
        by(1),-n(1),0,by(2),-n(2),0,by(3),-n(3),0,by(4),-n(4),0];
    k2 = k2 + 5/6*g*h*detJ*w* b'*b;
end
% Arrange stiffness in the usual u, v, w, ?x, ?y, ?z order for each node
nlm = [1, 2, 7, 8, 13, 14, 19, 20];
ke(nlm, nlm) = k1;
nlm = [3, 4, 5, 9, 10, 11, 15, 16, 17, 21, 22, 23];
ke(nlm, nlm) = k2;
rotzStiff = rotzStiffFactor*e*h*elemArea;
for i=6:6:24
    ke(i, i) = rotzStiff;
end

% compute the lumped stiffness matrix
masse = 0.25*elemArea*h;
rotme=.08333333*masse*h^2;
me = zeros(24,24);
me(1:12,1:12)=masse*eye(12);
me(13:24,13:24)=rotme*eye(12);

% Transform to global
ke = T'*ke*T;
re = T'*re;
end

