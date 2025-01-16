function N=shape_tetra10(xi)

% function N = shape_tetra10(xi)
%
%
%
r = xi(1); s = xi(2); t = xi(3);
rr = r^2;
ss = s^2;
tt = t^2;
rs = r*s;
st = s*t;
rt = t*r;
N(1) = 1 - 3*r + 2*rr + 4*rs + 4*rt - 3*s + 2*ss + 4*st - 3*t + 2*tt;
N(2) = -r + 2*rr;
N(3) = -s + 2*ss;
N(4) = -t + 2*tt;
N(5) = 4*r - 4*rr - 4*rs - 4*rt;
N(6) = 4*rs;
N(7) = -4*rs + 4*s - 4*ss - 4*st;
N(8) = -4*rt - 4*st + 4*t - 4*tt;
N(9) = 4*rt;
N(10)= 4*st;
