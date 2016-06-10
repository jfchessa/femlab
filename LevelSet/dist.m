function [d,xi,indx]=dist(p,p1,p2)

% *********************************************************************
% DIST: calculates the closest point from a point to a segment.  If the
%       normal distance is of the segment the endpoint is returned.

p=ones(size(p1,1),1)*p;

u=p2-p1;
v=p-p1;

dp=u(:,1).*v(:,1)+u(:,2).*v(:,2);
nv=sqrt(v(:,1).*v(:,1)+v(:,2).*v(:,2));
nu=sqrt(u(:,1).*u(:,1)+u(:,2).*u(:,2));

theta=zeros(size(nu));
nztheta=intersect(find(nu~=0),find(nv~=0));
theta(nztheta)=acos( dp(nztheta)./nu(nztheta)./nv(nztheta) );

dn=nv.*sin(theta);
tmp=p-p1;
d1=sqrt(tmp(:,1).*tmp(:,1)+tmp(:,2).*tmp(:,2));
tmp=p-p2;
d2=sqrt(tmp(:,1).*tmp(:,1)+tmp(:,2).*tmp(:,2));

xi=2*dp./nu-1;
d=dn;

caseIs=find( theta > pi/2 );
xi(caseIs)=-1;
d(caseIs)=d1(caseIs);

caseIIs=find( dp > nu.^2 );
xi(caseIIs)=1;
d(caseIIs)=d2(caseIIs);

[d,indx]=min(d);
xi=xi(indx,:);



