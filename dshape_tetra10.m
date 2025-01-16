function dNxi = dshape_tetra10(xi)

% function dNdxi = dshape_tetra10(xi)
%
%
%
r = xi(1); s = xi(2); t = xi(3);

dNxi(1,1) = -3 + 4*r + 4*s + 4*t;
dNxi(1,2) = dNxi(1,1);
dNxi(1,3) = dNxi(1,1);

dNxi(2,1) = -1 + 4*r;
dNxi(2,2) = 0.0;
dNxi(2,3) = 0.0;

dNxi(3,1) = 0.0;
dNxi(3,2) = -1 + 4*s;
dNxi(3,3) = 0.0;

dNxi(4,1) = 0.0;
dNxi(4,2) = 0.0;
dNxi(4,3) = -1 + 4*t;

dNxi(5,1) = 4 - 8*r - 4*s - 4*t;
dNxi(5,2) = -4*r;
dNxi(5,3) = -4*r;

dNxi(6,1) = 4*s;
dNxi(6,2) = 4*r;
dNxi(6,3) = 0.0;

dNxi(7,1) = -4*s;
dNxi(7,2) = 4 - 4*r - 8*s - 4*t;
dNxi(7,3) =-4*s;

dNxi(8,1) = -4*t;
dNxi(8,2) = -4*t;
dNxi(8,3) = 4 - 4*r - 4*s - 8*t;

dNxi(9,1) = 4*t;
dNxi(9,2) = 0.0;
dNxi(9,3) = 4*r;

dNxi(10,1) = 0.0;
dNxi(10,2) = 4*t;
dNxi(10,3) = 4*s;
