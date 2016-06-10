% tension.m
%
%
clear

node=[0.0 0.0;
      0.5 0.0;
      1.0 0.0;
      0.0 0.5;
      0.5 0.5;
      1.0 0.5;
      0.0 1.0;
      0.5 1.0;
      1.0 1.0;
      0.0 1.5;
      0.5 1.5;
      1.0 1.5];
    
element= [ 1  2  4;
           5  4  2;
           2  3  5;
           6  5  3;
           4  5  7;
           8  7  5;
           5  6  8;
           9  8  6;
           7  8 10;
          11 10  8;
           8  9 11;
          12 11  9 ];
        

young = 10E6;
poisson = 0.33;
thk=.125;

fext=zeros(24,1); 

ifix=[ 1 2 4 6 7 13 19 ]; % fix bottom edge in y and left edge in x
                                 
fext([22])=1000;   % point load on top in the y-direction                
fext([20 24])=500;



fea2d_tria3