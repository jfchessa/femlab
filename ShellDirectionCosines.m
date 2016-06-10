function H = ShellDirectionCosines(nodes)
% H = ShellDirectionCosines(nodes)
n1 = nodes(1,:);
n2 = nodes(2,:);
n3 = nodes(3,:);
L = sqrt((n2 - n1)*(n2 - n1)');
ex = (n2 - n1)/L;
ezz = cross(n2 - n1, n3 - n1);
Lt = sqrt(ezz*ezz');
ez = ezz/Lt;
ey = cross(ez, ex); 
H = [ex; ey; ez];
end