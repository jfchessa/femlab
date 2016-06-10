function  X=extrude_nodes(x1,y1,z1,v,nv)

% function X=extrude_nodes(x1,y1,z1,x2,y2)
%
%  extrudes the nodes given in x1, y1, z1, along the direction v with
%  nv nodes along the edge.  x1, y1, z1 must all be of the same length
%


nn=min([length(x1),length(y1),length(z1)]);

X=zeros(nn*nv,3);

dv=v/(nv-1);
n=1;
for r=0:nv-1
    
    X(n:n+nn-1,1)=x1+r*dv(1);
    X(n:n+nn-1,2)=y1+r*dv(2);
    X(n:n+nn-1,3)=z1+r*dv(3);
    
    n=n+nn;
    
end
