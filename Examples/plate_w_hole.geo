lc = 0.01;
lc2=.5;
W = 5;
H=10;
R=.25;

Point(1)={0.0,0.0,0.0,lc};
Point(2)={ R,0.0,0.0,lc};
Point(3)={R*Sqrt(2)/2,R*Sqrt(2)/2,0.0,lc};
Point(4)={0.0,R,0.0,lc};

Point(5)={W/2,0.0,0.0,lc2};
Point(6)={W/2,H/2,0.0,lc2};
Point(7)={0.0,H/2,0.0,lc2};
Line(1) = {2,5};
Line(2) = {5,6};
Line(3) = {6,7};
Line(4) = {7,4};
Circle(5) = {2,1,3};
Circle(6) = {3,1,4};


Line Loop(7) = {2,3,4,-6,-5,1};
Plane Surface(8) = {7};

Physical Line(9) = {1};    // bottom edge
Physical Line(10) = {4};   // left edge
Physical Line(11) = {3};   // top edge
Physical Surface(12) = {8};  // domain
