clear

coord=[0 0 0;0 0 5;0 0 1];
[pts,wts]=quadrature('line3',5);

l=0;
for i=1:length(pts)
    
    disp('----------------------------')
    dNdxi=dshape_line3(pts(i,:))
    [jac,jmat]=element_jacobian(coord,dNdxi)
    l=l+jac*wts(i)
    
end

disp('----------------------------')
fprintf('the length of the line element is %5.2f\n', l);

norm( coord(2,:)-coord(1,:) )

coord=[0 0 0; 1 0 0; 0 1 0.5];
[pts,wts]=quadrature('tria3',5);
A=0;
for i=1:length(pts)
    
    disp('----------------------------')
    dNdxi=dshape_tria3(pts(i,:))
    [jac,jmat]=element_jacobian(coord,dNdxi)
    A=A+jac*wts(i)
    
end

disp('----------------------------')
fprintf('the area of the traingle element is %5.2f\n', A);

