function ke=kmat_truss2d(coord,AE)

% function ke=kmat_truss2d(coord,AE)
%
% Computes the element stiffness matrix for a 2D linear elastic
% truss element.  
%
%  coord - a 2 by 2 matrix of the node coordinates for the element
%  AE    - the product of the cross sectional area times the Young's
%           modulus for the element
%
% Jack Chessa jfchessa@utep.edu
%

le = norm( coord(2,:)-coord(1,:) );   % length of the element

c=(coord(2,1)-coord(1,1))/le;  % element cosine
s=(coord(2,2)-coord(1,2))/le;  % element sine

c2 = c^2; % cos^2
s2 = s^2; % sin ^2
cs = c*s; % cos*sin 

ke= AE/le*[  c2  cs -c2 -cs;
             cs  s2 -cs -s2;
            -c2 -cs  c2  cs;
            -cs -s2  cs  s2 ];
          