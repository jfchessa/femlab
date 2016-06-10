function [B,A]=bmat_quad8(coord,xi)

% function B=bmat_quad8(coord,xi)
%
% Computes the strain-displacement matrix (B matrix) for a eight node
% quadrilateral element.
%
%    coord: the nodal coordinates of the element (4x2 matrix)
%
% function [B,jac]=bmat_quad4(coord,xi)
%
% Computes the B matrix and the Jacobian
%
% Written by Jack Chessa, jfchessa@utep.edu


r=0.5*(xi(1)+1.0)
s=0.5*(xi(2)+1.0)

dNdxi = [ ( s - 1.0  ) * ( - 4.0  * r - 2.0  * s + 3.0  ),...
  ( r - 1.0  ) * ( - 4.0  * s - 2.0  * r + 3.0  );
  ( s - 1.0  ) * ( - 4.0  * r + 2.0  * s + 1.0  ),...
  r *       (   4.0  * s - 2.0  * r - 1.0  );
  s         * (   4.0  * r + 2.0  * s - 3.0  ),...
  r *       (   4.0  * s + 2.0  * r - 3.0  );
  s         * (   4.0  * r - 2.0  * s - 1.0  ),...
  ( r - 1.0  ) * ( - 4.0  * s + 2.0  * r + 1.0  );
  4.0  * ( 2.0  * r - 1.0  )     * ( s - 1.0  ),...
  4.0  * r * ( r - 1.0  );
  - 4.0  *                     s * ( s - 1.0  ),...
  - 4.0  * r               * ( 2.0  * s - 1.0  );
  - 4.0  * ( 2.0  * r - 1.0  ) * s,...
  - 4.0  * r * ( r - 1.0  );
  4.0  *                     s * ( s - 1.0  ),...
  4.0  *     ( r - 1.0  ) * ( 2.0  * s - 1.0  ) ];

J=coord'*dNxi;
jac=det(J);
dN=dNxi*inv(J);

B=zeros(3,16);
for i=1:8
  B(1,2*i-1)=dN(i,1);
  B(2,2*i)=dN(i,2);
  B(3,2*i-1)=dN(i,2);
  B(3,2*i)=dN(i,1);
end