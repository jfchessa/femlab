function [ke, rq] = kmat_beam3d(coord, ptoffa, e, G, Ir, Is, J, A, qr, qs )

% [ke, rq] = kmat_beam3d(coord, E, G, Ir, Is, J, A, qr, qs )
% 	Generates equations for a space frame element in 3D
%
% 	coord = coordinates at the element ends
% 	ptoffa =  point of of beam axis to define qr plane
% 	E = modulus of elasticity
% 	G = shear modulus
% 	Ir, Is = moment of inertias about element r and s axes
% 	J = torsional rigity
% 	A = area of cross-section
% 	qr, qs = distributed loads along the element r and s axes (optional)

EIr=e*Ir; EIs=e*Is; GJ=G*J; EA = e*A;
n1=coord(1,1:3); n2=coord(2,1:3); n3=ptoffa;
L=sqrt(dot((n2-n1),(n2-n1)));
ex = (n2 - n1)/L;
eyy = cross(n3 - n1, n2 - n1);
ey = eyy/sqrt(dot(eyy,eyy));
ez = cross(ex, ey);
H = [ex; ey; ez];
T = zeros(12);
T([1, 2, 3], [1, 2, 3]) = H; 
T([4,5,6], [4,5,6]) = H;
T([7,8,9], [7,8,9]) = H;
T([10,11,12], [10,11,12]) = H;
TT = T';
ke = [EA/L, 0, 0, 0, 0, 0, -(EA/L), 0, 0, 0, 0, 0;
    0, (12*EIr)/L^3, 0, 0, 0, (6*EIr)/L^2, 0, -((12*EIr)/L^3), ...
        0,0, 0, (6*EIr)/L^2;
    0, 0, (12*EIs)/L^3, 0, -((6*EIs)/L^2), 0, ...
        0, 0, -((12*EIs)/L^3), 0, -((6*EIs)/L^2), 0;
    0, 0, 0, GJ/L, 0, 0, 0, 0, 0, -(GJ/L), 0, 0;
    0, 0, -((6*EIs)/L^2), 0, (4*EIs)/L, 0, 0, 0, (6*EIs)/L^2, ...
        0,(2*EIs)/L, 0;
    0, (6*EIr)/L^2, 0, 0, 0, (4*EIr)/L, 0, ...
        -((6*EIr)/L^2), 0, 0, 0, (2*EIr)/L;
    -(EA/L), 0, 0, 0, 0, 0, EA/L, 0, 0, 0, 0, 0;
    0, -((12*EIr)/L^3), 0, 0, 0, -((6*EIr)/L^2), ...
        0, (12*EIr)/L^3, 0, 0, 0, -((6*EIr)/L^2);
    0, 0, -((12*EIs)/L^3), 0, (6*EIs)/L^2, 0, 0, 0, ...
        (12*EIs)/L^3, 0, (6*EIs)/L^2, 0;
    0, 0, 0, -(GJ/L), 0, 0, 0, 0, 0, GJ/L, 0, 0;
    0, 0, -((6*EIs)/L^2), 0, (2*EIs)/L, 0, 0,0, (6*EIs)/L^2, ...
        0, (4*EIs)/L, 0;
    0, (6*EIr)/L^2, 0, 0, 0, (2*EIr)/L, 0, -((6*EIr)/L^2), 0, ...
        0, 0, (4*EIr)/L];
ke = TT*ke*T;
if ( nargin > 8 )
	rq = TT*[0; (L*qs)/2; (L*qr)/2; 0; -((L^2*qr)/12); ...
        (L^2*qs)/12; 0; (L*qs)/2;
    	(L*qr)/2; 0; (L^2*qr)/12; -((L^2*qs)/12)];
else
	rq=zeros(12,1);
end


