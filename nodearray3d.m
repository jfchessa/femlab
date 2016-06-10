function X=nodearray3d( corners, n1, n2, n3, st1, sp1, st2, sp2, st3, sp3 )

% function X=NODEARRAY3D( CORNERS, N1, N2, N3 )
%
% Gerenates an array of nodes on a 3D domain
%
% function X=NODEARRAY3D( CORNERS, N1, N2, N3, SCALETYPE1, SCALEPARAM1,
%             SCALETYPE2, SCALEPARAM2, SCALETYPE3, SCALEPARAM3 )
%

if ( nargin <5 )
    st1='NONE';
    sp1=1;
    st2='NONE';
    sp2=1;
    st3='NONE';
    sp3=1;
end


if ( size(corners,1)==8 )
    shapefun='shape_hexa8';
elseif ( size(corners,1)==4 )
    shapefun='shape_tetra4';
else 
    error('incorrect dimension for CORNERS in NODEARRAY3D');
end

xi=nodearray1d(-1,1,n1,st1,sp1);
eta=nodearray1d(-1,1,n2,st2,sp2);
zeta=nodearray1d(-1,1,n3,st3,sp3);

ni=length(xi); nj=length(eta); nk=length(zeta);
nn=ni*nj*nk;
nd=size(corners,2);

X=zeros(nn,3);

n=1;
for k=1:nk
    for j=1:nj
        for i=1:ni
            pt=[xi(i),eta(j),zeta(k)];
            N=shape_hexa8(pt)';
            X(n,1:nd)=N*corners;
            n=n+1;
        end
    end
end
