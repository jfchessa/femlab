% frame2d_ex.m
%
% A simple 2D frame example
%

clear

node=[ 0 0; 0 10; 4 10; 4 0 ];
conn=[ 1 2; 2 3; 3 4 ];

A=1.0;
E=10e6;
G=3E6;
I=2;

nn=size(node,1);
ne=size(conn,1);
ndof=nn*3;
K=sparse(ndof,ndof);

ifix=[1:3 10:12];
f=zeros(ndof,1);
f(4)=10;

for e=1:size(conn,1)
    
    sctr = get_scatter( conn(e,:), 3 );
    K(sctr,sctr)=K(sctr,sctr)+kmat_beam( node(conn(e,:),:), A, E, I );
    
end
    
[d,freac]=fesolve(K,f,ifix);