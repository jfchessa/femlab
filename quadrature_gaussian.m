function [ quadpoint, quadweight ] = quadrature_gaussian(npts,sdim)

% function [ QPTS, QWTS ] = QUADRATURE_GAUSSIAN(NPTS,DIM)
%
% Returns the quadrature points and weights for gaussian quadrature
%
% NPTS - The polynomial order of the integrand
% DIM - Spacial dimension to be integrated (default is 1)
%
% This is part of FEMLAB
% Written by Jack Chessa, jfchessa@utep.edu

if ( nargin==1 )
    sdim=1;
end
% 
% if ( quadorder > 15 )  % check for valid quadrature order
%     disp('Order of quadrature to high for QUADRATURE_GAUSSIAN');
%     quadorder = 15;
% end
% 
% npts=ceil( 0.5*(quadorder+1) );

quadpoint=zeros(npts,1); quadweight=zeros(npts,1);

switch ( npts )
    case 1
        quadpoint(1) = 0.000000000000000;
        quadweight(1) = 2.000000000000000;
        
    case 2
        quadpoint(1) = 0.577350269189626;
        quadpoint(2) =-0.577350269189626;
        
        quadweight(1) = 1.000000000000000;
        quadweight(2) = 1.000000000000000;
        
    case 3
        quadpoint(1) = 0.774596669241483;
        quadpoint(2) =-0.774596669241483;
        quadpoint(3) = 0.000000000000000;
        
        quadweight(1) = 0.555555555555556;
        quadweight(2) = 0.555555555555556;
        quadweight(3) = 0.888888888888889;
        
    case 4
        quadpoint(1) = 0.861134311594053;
        quadpoint(2) =-0.861134311594053;
        quadpoint(3) = 0.339981043584856;
        quadpoint(4) =-0.339981043584856;
        
        quadweight(1) = 0.347854845137454;
        quadweight(2) = 0.347854845137454;
        quadweight(3) = 0.652145154862546;
        quadweight(4) = 0.652145154862546;
        
    case 5
        quadpoint(1) = 0.906179845938664;
        quadpoint(2) =-0.906179845938664;
        quadpoint(3) = 0.538469310105683;
        quadpoint(4) =-0.538469310105683;
        quadpoint(5) = 0.000000000000000;
        
        quadweight(1) = 0.236926885056189;
        quadweight(2) = 0.236926885056189;
        quadweight(3) = 0.478628670499366;
        quadweight(4) = 0.478628670499366;
        quadweight(5) = 0.568888888888889;
        
    case 6
        quadpoint(1) = 0.932469514203152;
        quadpoint(2) =-0.932469514203152;
        quadpoint(3) = 0.661209386466265;
        quadpoint(4) =-0.661209386466265;
        quadpoint(5) = 0.238619186003152;
        quadpoint(6) =-0.238619186003152;
        
        quadweight(1) = 0.171324492379170;
        quadweight(2) = 0.171324492379170;
        quadweight(3) = 0.360761573048139;
        quadweight(4) = 0.360761573048139;
        quadweight(5) = 0.467913934572691;
        quadweight(6) = 0.467913934572691;
        
    case 7
        quadpoint(1) =  0.949107912342759;
        quadpoint(2) = -0.949107912342759;
        quadpoint(3) =  0.741531185599394;
        quadpoint(4) = -0.741531185599394;
        quadpoint(5) =  0.405845151377397;
        quadpoint(6) = -0.405845151377397;
        quadpoint(7) =  0.000000000000000;
        
        quadweight(1) = 0.129484966168870;
        quadweight(2) = 0.129484966168870;
        quadweight(3) = 0.279705391489277;
        quadweight(4) = 0.279705391489277;
        quadweight(5) = 0.381830050505119;
        quadweight(6) = 0.381830050505119;
        quadweight(7) = 0.417959183673469;
        
    case 8
        quadpoint(1) =  0.960289856497536;
        quadpoint(2) = -0.960289856497536;
        quadpoint(3) =  0.796666477413627;
        quadpoint(4) = -0.796666477413627;
        quadpoint(5) =  0.525532409916329;
        quadpoint(6) = -0.525532409916329;
        quadpoint(7) =  0.183434642495650;
        quadpoint(8) = -0.183434642495650;
        
        quadweight(1) = 0.101228536290376;
        quadweight(2) = 0.101228536290376;
        quadweight(3) = 0.222381034453374;
        quadweight(4) = 0.222381034453374;
        quadweight(5) = 0.313706645877887;
        quadweight(6) = 0.313706645877887;
        quadweight(7) = 0.362683783378362;
        quadweight(8) = 0.362683783378362;
        
    otherwise
        disp('Order of quadrature to high for QUADRATURE_GAUSSIAN');
        
end  % end of quadorder switch


if ( sdim == 2 )
    [quadpoint,quadweight] = quadrature_compound(quadpoint,quadweight,...
            quadpoint,quadweight);
    
elseif ( sdim == 3 )
    [quadpoint,quadweight] = quadrature_compound(quadpoint,quadweight,...
            quadpoint,quadweight,quadpoint,quadweight);

    
end



end

