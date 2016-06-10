function [bmat,rpt,jac,N] = bmat_axisymmetric( coord, xi, type )

% function [bmat,r,A,N] = bmat_axisymmetric( coord, xi, type )
%
%   Computes the B-matrix for an axisymmetric element. For the axisymmetric
%   case the order of the strains are as follows:
%           [ rr, theta, zz, rz ]
%
%
%       ^  z
%       |
%       |      
%       |
%       |    x (r,z)
%       |
%       |
%       |________________>  r
%
%   This element assumes the x axis is the radial axis and the axis of 
%   symmetry is the z axis.
%
%   coord - the nodal coordinates of the element (nnx2 matrix)
%   xi    - the coordinate in the parent element space 
%   type  - 'quad4', 'tria3', .. or any valid 2d element (calls the 
%           associated shape_ and bmat_ functions.
%
%   bmat - is the axisymmetric B matrix
%     r - is the radius at the parent coord xi.
%   jac -  is the determinate of the element (or area for simplex elements)
%   N -  is the element shape function vector
%
% -------------------
%
%     Written by Jack Chessa, jfchessa@utep.edu
%     for the FEMLAB library
%

N = feval( ['shape_',type], xi );
rpt = N'*coord(:,1);

nn=length(N);
bmat=zeros(4,2*nn);

[ bmat([1 3 4],:), jac] = feval(['bmat_',type], coord, xi );
bmat( 2, 1:2:2*nn ) = N'./rpt;



