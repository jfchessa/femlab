function [ qpts, qwts ] = quadrature_compound(p1,w1,p2,w2,p3,w3)

% function [ QPTS, QWTS ] = QUADRATURE_COMPOUND(P1,W1,P2,W2,P3,W3)
%
% Computes a compound quadrature rule.  This can be used generate 2D and
% 3dD quadrature rules from tensorial products of 1D and 2D quadrature
% rules.
%

n=1;
if ( nargin==4 )
    
    npts=length(w1)*length(w2);
    sdim=size(p1,2)+size(p2,2);
    qwts=zeros(npts,1);
    qpts=zeros(npts,sdim);
    for j = 1:length(w2)
        for i = 1:length(w1)
            qpts(n,:) = [ p1(i,:), p2(j,:)];
            qwts(n) = w1(i)*w2(j);
            n = n+1;
        end
    end
    
else % 3 quadrature rules
    
    npts=length(w1)*length(w2)*length(w3);
    sdim=size(p1,2)+size(p2,2)+size(p3,2);
    qwts=zeros(npts,1);
    qpts=zeros(npts,sdim);
    for k = 1:length(w3)
        for j = 1:length(w2)
            for i = 1:length(w1)
            qpts(n,:) = [ p1(i,:), p2(j,:) p3(k,:) ];
            qwts(n) = w1(i)*w2(j)*w3(k);
            n = n+1;
        end
    end    
    
end


end
