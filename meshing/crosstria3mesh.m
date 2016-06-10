function [node,conn]=crosstria3mesh(corners,he,st1,sp1,st2,sp2,plotit)


if ( nargin < 3 ), st1='NONE'; end
if ( nargin < 4 ), sp1=1; end
if ( nargin < 5 ), st2='NONE'; end
if ( nargin < 6 ), sp2=1; end
if (nargin<7), plotit=0;  end

length = 0.5*( norm(corners(1,:)-corners(2,:)) + ...
      norm(corners(3,:)-corners(4,:)) );
height = 0.5*( norm(corners(3,:)-corners(2,:)) + ...
      norm(corners(1,:)-corners(4,:)) );

nnx=2*ceil(2*length/he)+1;
nny = 2*ceil(2*height/he)+1;

node=nodearray2d( corners, nnx, nny, st1, sp1, st2, sp2 );
                

conn=genconn2d([nnx+1 1 nnx+2],(nnx-1)/2,(nny-1)/2,2,3+nnx);
conn=[conn; genconn2d([1 2 nnx+2],(nnx-1)/2,(nny-1)/2,2,3+nnx)];
conn=[conn; genconn2d([2 3 nnx+2],(nnx-1)/2,(nny-1)/2,2,3+nnx)];
conn=[conn; genconn2d([3 nnx+3 nnx+2],(nnx-1)/2,(nny-1)/2,2,3+nnx)];
conn=[conn; genconn2d([nnx+3 2*nnx+3 nnx+2],(nnx-1)/2,(nny-1)/2,2,3+nnx)];
conn=[conn; genconn2d([2*nnx+3 2*nnx+2 nnx+2],(nnx-1)/2,(nny-1)/2,2,3+nnx)];
conn=[conn; genconn2d([2*nnx+2 2*nnx+1 nnx+2],(nnx-1)/2,(nny-1)/2,2,3+nnx)];
conn=[conn; genconn2d([2*nnx+1 nnx+1 nnx+2],(nnx-1)/2,(nny-1)/2,2,3+nnx)];

if ( plotit )
    trisurf(conn,node(:,1),node(:,2),0*node(:,1))
    view(2)
    axis equal
end
