function X=nodearray2d( corners, n1, n2, st1, sp1, st2, sp2 )

% function X=NODEARRAY2D( CORNERS, N1, N2 )
%
% Gerenates an array of nodes on a 2D region which may be in 3D
%
% function X=NODEARRAY2D( CORNERS, N1, N2, SCALETYPE1, SCALEPARAM1,
%             SCALETYPE2, SCALEPARAM2 )
% CORNERS is an array of the coordinates of the four corner points.  Each
% row is a point and the order of teh points should be in a right handed
% sence.

if ( nargin <4 )
    st1='NONE';
    sp1=1;
    st2='NONE';
    sp2=1;
end


if ( size(corners,1)==4 )
    shapefun='shape_quad4';
elseif ( size(corners,1)==3 )
    shapefun='shape_tria3';
else 
    error('incorrect dimension for CORNERS in NODEARRAY2D');
end

xi=nodearray1d(-1,1,n1,st1,sp1);
eta=nodearray1d(-1,1,n2,st2,sp2);

ni=length(xi); nj=length(eta); 
nn=ni*nj;
nd=size(corners,2);

X=zeros(nn,nd);

n=1;
for j=1:nj
    for i=1:ni
        pt=[xi(i),eta(j)];
        N=shape_quad4(pt)';
        X(n,1:nd)=N*corners;
        n=n+1;
    end
end
