%function [bmat,N] = bmat_axisymmetric( xi, ri, type )

% function 
%
%
% type: 'quad4', 'tria3'
%
%
%  Note: for the axisymmetriccase the order of the strains are as follows
%       [ rr, theta, zz, rz ]

clear;
coord=[1/3,1/3];
ri=[0;1;0];

dNxi = dshape_tria3(coord);
N  = shape_tria3(coord);
r  = N'*ri;

nn=length(N);
bmat=zeros(4*nn,2);

for i=1:nn
    
    coli=(2*i-1):2*i;
    
    bmat( :, coli ) = [ dN(i,1),     0.0; 
                         N(i)/r,     0.0; 
                            0.0, dN(i,2); 
                        dN(i,2), dN(i,1) ];   
    
end




