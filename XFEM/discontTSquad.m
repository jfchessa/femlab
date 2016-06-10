function [W,Q]=discontTSquad(order,phi1s,phi2s,plot_quad)

% function discontTSquad
%
% [W,Q]=discont_quad(ORDER,PHI)
%
%      equad returns the quadratrue weights and points for integration of
%      discontinous triangular elements in the parent coordinate system of
%      the element.  The location of the lines of discontinuity are defined 
%      by PHI a vector of nodal signed distances to the discontinuity. ORDER
%      is the polynomial order for wich the quadrature will be exact.
%
%  [W,Q]=discont_quad(ORDER,PHI1,PHI2)
%   
%      allows for the possibility of two lines of discontiuity in a 
%      single element as given by PHI1 and PHI2 respectivly.
%
%      This function splits the element in to subelements that conform
%      to the discontinuities, and  map the appropriate quadrature rules
%      from the subelements onto the parent coordinate space of the 
%      initial element.
%
%   [W,Q]=discontT3quad(order,phi1s,phi2s,1)
%      Will perform the quadratre and plot the subelements and locations of 
%      the quadrature points (cool)


% clear
% order=1;
% num=2;
% phi1s=[.1;-.4;.1];
% phi2s=[.2; -.3; .2];
% plot_quad=1;
% END OF DEGUG VALUES

if ( nargin == 2 )
  phi2s=phi1s;
  num=1;
else
  num=2;
end

if( nargin < 4 )
  plot_quad=0;
end

% check if phis do not cut triangle
if ( max(phi1s)*min(phi1s) >= 0 & max(phi2s)*min(phi2s) >= 0  ) 
  [W,Q]=quadrature(order,'TRIANGULAR',2);
  return
else
  
  % split triangle with phis
  edge=[1 2 3 1];
  edgeXi=[0 0;1 0;0 1;0 0];
  normal=[0 -1;1 1;-1 0];
  plump=.0000001;
  nodes=edgeXi(1:3,:);
  
  phis=phi1s;
  for p=1:num
    
    for e=1:3  % loop over triangle edges
      n1=edge(e);
      n2=edge(e+1);
      
      if ( phis(n1)*phis(n2)<0 ) % edge is split
        xi=interp1( phis(edge(e:e+1)), edgeXi(e:e+1,1), 0 );
        eta=interp1( phis(edge(e:e+1)), edgeXi(e:e+1,2), 0 );
        nodes=[nodes;[xi,eta]+plump*normal(e)];
      end
      
    end
    phis=phi2s;
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if ( size(nodes,1) == 7 )
    
    nadd=4;
    addedNodes1 = [linspace(nodes(4,1),nodes(5,1),nadd);
                   linspace(nodes(4,2),nodes(5,2),nadd)]';
    
    addedNodes2 = [linspace(nodes(6,1),nodes(7,1),nadd);
                   linspace(nodes(6,2),nodes(7,2),nadd)]';
    
    nodes=[nodes;addedNodes1(2:nadd-1,:);addedNodes2(2:nadd-1,:)];
    
  end
  
  % get decompused triangles
  tri=delaunay(nodes(:,1),nodes(:,2));
  tri=tricheck(nodes,tri);
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % loop over subtriangles to get quadrature points and weights
  pt=1;
  for e=1:size(tri,1)
    
    [w,q]=quadrature(order,'TRIANGULAR',2);
    % transform quadrature points into the parent element
    coord=nodes(tri(e,:),:);
    a=det([coord,[1;1;1]])/2;
    if ( a<0 )  % need to swap connectivity
      coord=[coord(2,:);coord(1,:);coord(3,:)];
      a=det([coord,[1;1;1]])/2;
    end
    
    if ( a~=0 )
      for n=1:length(w)
        N=lagrange_basis('T3',q(n,:));
        Q(pt,:)=N'*coord;
        W(pt,1)=2*w(n)*a;
        pt=pt+1;
      end
    end
    
  end
end

if ( plot_quad )
  clf; 
  trisurf(tri,nodes(:,1),nodes(:,2),0*nodes(:,1)); 
  view(2)
  hold on
  plot( Q(:,1), Q(:,2), 'rx' )
end


    
    
    
    