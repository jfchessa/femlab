function B=fill_bmatrix(dNx,sdim,nn)

% function B=fill_bmatrix(dNx,sdim,nn)
%
% Given the gradient of the shape functions w.r.t. the spacial coordinates
% fills out a B-matrix.
%

if ( nargin<2 )
    sdim=size(dNx,2);
end
if ( nargin<3 )
    nn=size(dNx,1);
end

if ( sdim==1 )
    
    B=dNx(1:nn,1)';  
    
elseif ( sdim==2 )
    
    B=zeros(3,sdim*nn);
    B(1,1:2:sdim*nn)=dNx(1:nn,1)'; 
    B(2,2:2:sdim*nn)=dNx(1:nn,2)'; 
    B(3,1:2:sdim*nn)=B(2,2:2:sdim*nn); 
    B(3,2:2:sdim*nn)=B(1,1:2:sdim*nn);
    
else
    
    B=zeros(6,sdim*nn);
    B(1,1:3:sdim*nn)=dNx(1:nn,1)'; 
    B(2,2:3:sdim*nn)=dNx(1:nn,2)';  
    B(3,3:3:sdim*nn)=dNx(1:nn,3)'; 
    B(4,2:3:sdim*nn)=B(3,3:3:sdim*nn);
    B(4,3:3:sdim*nn)=B(2,2:3:sdim*nn);
    B(5,3:3:sdim*nn)=B(1,1:3:sdim*nn);
    B(5,1:3:sdim*nn)=B(3,3:3:sdim*nn);
    B(6,1:3:sdim*nn)=B(2,2:3:sdim*nn);
    B(6,2:3:sdim*nn)=B(1,1:3:sdim*nn);
    
end