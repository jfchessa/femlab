function [ quadpoint, quadweight ] = quadrature_simplex(quadorder,sdim)

% function [ QPTS, QWTS ] = QUADRATURE_SIMPLEX(ORDER,DIM)
%
% Returns the quadrature points and weights for simplex elements (ie.
% triangular and tetrahedral elements)
%
% ORDER - The polynomial order of the integrand
% DIM - Spacial dimension to be integrated (default is 2)
%
% This is part of FEMLAB
% Written by Jack Chessa, jfchessa@utep.edu

if ( nargin==1 )
    sdim=2;
end

if ( sdim == 3 )  %%% TETRAHEDRA
      
      if ( quadorder ~= 1 &  quadorder ~= 2 &  quadorder ~= 3  ) 
        % check for valid quadrature order
        disp('Incorect quadrature order for QUADRATURE_SIMPLEX');
        quadorder = 1;
      end
      
      if  ( quadorder == 1 )
        quadpoint = [ 0.25 0.25 0.25 ];
        quadweight = 1;
        
      elseif ( quadorder == 2 ) 
        quadpoint = [ 0.58541020  0.13819660  0.13819660;
                      0.13819660  0.58541020  0.13819660;
                      0.13819660  0.13819660  0.58541020;
                      0.13819660  0.13819660  0.13819660];
        quadweight = [1; 1; 1; 1]/4;
        
      elseif ( quadorder == 3 ) 
        quadpoint = [ 0.25  0.25  0.25;
                      1/2   1/6   1/6;
                      1/6   1/2   1/6;
                      1/6   1/6   1/2;
                      1/6   1/6   1/6];
        quadweight = [-4/5 9/20 9/20 9/20 9/20]';
        
      end
      
      quadweight=0.166666666666666666666667*quadweight;
         
    else  %%% TRIANGLES
      
      if ( quadorder > 7 ) % check for valid quadrature order
        disp('Quadrature order too high for QUADRATURE_SIMPLEX');
        quadorder = 1;
      end
      
      if ( quadorder <= 1 )   % set quad points and quadweights
        quadpoint = [ 0.3333333333333, 0.3333333333333 ];
        quadweight = 1;
        
      elseif ( quadorder == 2 ) 
        quadpoint = zeros( 3, 2 );
        quadweight = zeros( 3, 1 );
        
        quadpoint(1,:) = [ 0.1666666666667, 0.1666666666667 ];
        quadpoint(2,:) = [ 0.6666666666667, 0.1666666666667 ];
        quadpoint(3,:) = [ 0.1666666666667, 0.6666666666667 ]; 
        
        quadweight(1) = 0.3333333333333; 
        quadweight(2) = 0.3333333333333; 
        quadweight(3) = 0.3333333333333;   
        
      elseif ( quadorder <= 5 ) 
        quadpoint = zeros( 7, 2 );
        quadweight = zeros( 7, 1 );
        
        quadpoint(1,:) = [ 0.1012865073235, 0.1012865073235 ];
        quadpoint(2,:) = [ 0.7974269853531, 0.1012865073235 ];
        quadpoint(3,:) = [ 0.1012865073235, 0.7974269853531 ]; 
        quadpoint(4,:) = [ 0.4701420641051, 0.0597158717898 ];
        quadpoint(5,:) = [ 0.4701420641051, 0.4701420641051 ];
        quadpoint(6,:) = [ 0.0597158717898, 0.4701420641051 ]; 
        quadpoint(7,:) = [ 0.3333333333333, 0.3333333333333 ];
        
        quadweight(1) = 0.1259391805448; 
        quadweight(2) = 0.1259391805448; 
        quadweight(3) = 0.1259391805448; 
        quadweight(4) = 0.1323941527885;
        quadweight(5) = 0.1323941527885;
        quadweight(6) = 0.1323941527885;
        quadweight(7) = 0.2250000000000;  
        
      else
        quadpoint = zeros( 13, 2 );
        quadweight = zeros( 13, 1 );
        
        quadpoint(1 ,:) = [ 0.0651301029022, 0.0651301029022 ];
        quadpoint(2 ,:) = [ 0.8697397941956, 0.0651301029022 ];
        quadpoint(3 ,:) = [ 0.0651301029022, 0.8697397941956 ];
        quadpoint(4 ,:) = [ 0.3128654960049, 0.0486903154253 ];
        quadpoint(5 ,:) = [ 0.6384441885698, 0.3128654960049 ];
        quadpoint(6 ,:) = [ 0.0486903154253, 0.6384441885698 ];
        quadpoint(7 ,:) = [ 0.6384441885698, 0.0486903154253 ];
        quadpoint(8 ,:) = [ 0.3128654960049, 0.6384441885698 ];
        quadpoint(9 ,:) = [ 0.0486903154253, 0.3128654960049 ];
        quadpoint(10,:) = [ 0.2603459660790, 0.2603459660790 ];
        quadpoint(11,:) = [ 0.4793080678419, 0.2603459660790 ];
        quadpoint(12,:) = [ 0.2603459660790, 0.4793080678419 ];
        quadpoint(13,:) = [ 0.3333333333333, 0.3333333333333 ];
        
        quadweight(1 ) = 0.0533472356088;
        quadweight(2 ) = 0.0533472356088; 
        quadweight(3 ) = 0.0533472356088;
        quadweight(4 ) = 0.0771137608903;
        quadweight(5 ) = 0.0771137608903;
        quadweight(6 ) = 0.0771137608903;
        quadweight(7 ) = 0.0771137608903;
        quadweight(8 ) = 0.0771137608903;
        quadweight(9 ) = 0.0771137608903;
        quadweight(10) = 0.1756152576332; 
        quadweight(11) = 0.1756152576332; 
        quadweight(12) = 0.1756152576332;
        quadweight(13) =-0.1495700444677; 
        
      end 

      quadweight=0.5*quadweight;
    end


end

