function C = cmat_mat1(E, nu, form)

% function C = CMAT_MAT1(E, NU, FORMULATION)
% 
% Computes the material stiffness matrix for a linear isotropic elastic
% material (like a NASTRAN MAT1)
%
% E = Youngs modulus
% nu = Poissons ratio
% formulation = 'PSTRESS', 'PSTRAIN', '3D', 'AXISYMMETRIC', default is '3D'
%
%  Note: for the axisymmetric case the order of the strains are as follows
%       [ rr, theta, zz, rz ]

if ( nargin<3 )
    form='3D';
end

switch lower(form)

    case 'pstress'

        c1=E/(1-nu^2);   
        c2=nu*c1;
        c3=0.5*(1-nu)*c1;

        C=[c1 c2 0;c2 c1 0; 0 0 c3];
        
    case 'pstrain'
        
        c0=E/(1-2*nu)/(1+nu);
        c1=(1-nu)*c0;   
        c2=nu*c0;   
        c3=0.5*(1-2*nu)*c0;

        C=[c1 c2 0;c2 c1 0; 0 0 c3];
        
    case 'axisymmetric'
               
        c0=E/(1-2*nu)/(1+nu);
        c1=(1-nu)*c0;   
        c2=nu*c0;   
        c3=0.5*(1-2*nu)*c0;

        C=[ c1 c2 c2 0;
            c2 c1 c2 0; 
            c2 c2 c1 0; 
             0  0  0 c3];
        
    otherwise % '3d'

        c0=E/(1-2*nu)/(1+nu);
        c1=(1-nu)*c0;   
        c2=nu*c0;   
        c3=0.5*(1-2*nu)*c0;

        C=[c1 c2 c2 0 0 0;
           c2 c1 c2 0 0 0; 
           c2 c2 c1 0 0 0;
            0  0  0 c3 0 0; 
            0  0  0 0 c3 0; 
            0  0  0 0  0 c3 ];
        
end


