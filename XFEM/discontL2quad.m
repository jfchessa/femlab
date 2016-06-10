function [W,Q]=equadL2(order,phi1s,phi2s)

% enriched quadrature for L2 elements

% clear
% order=3;
% phi1s=[-1 3]';
% phi2s=phi1s;

if ( nargin == 2 )
  phi2s=phi1s;
  num=1;
else
  num=2;
end

if ( phi1s==phi2s )
  num=1;
end

if ( phi1s(1)*phi1s(2) >= 0 & phi2s(1)*phi2s(2) >= 0 )
  
  [W,Q]=quadrature(order,'GAUSS',1);
  
else
  
  % split element with zero points from phis
  phis=phi1s;
  nodes=[-1;1];
  W=[];
  Q=[];
  for n=1:num
    if ( phis(1)*phis(2) < 0 )
      xi=( phis(1) + phis(2) )/( phis(1) - phis(2) );
      nodes=[nodes;xi];
    end
    phis=phi2s;
  end
  nodes=sort(nodes);
  
  for n=1:(size(nodes,1)-1)
    [w,q]=quadrature(order,'GAUSS',1);
    l=nodes(n+1)-nodes(n);
    W=[W;w*l/2];
    Q=[Q; [(1-q)/2, (1+q)/2]*nodes([n,n+1]) ];
  end
    
end

%plot(Q,W,'ro')