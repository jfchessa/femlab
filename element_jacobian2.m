function [jac,jmat]=element_jacobian(dNdxi,coord,sdim,edim)

% function [J,JMAT]=ELEMENT_JACOBIAN(DNDXI,COORD,SDIM,EDIM)
%
% Computes the Jacobian, jac, and the Jacobian matrix, jmat, for an
% isoparametric element.  This is the Jacobian of the transformation
% betweent the parent coordinate and spacial coordinate systems.
%
%    DNDXI -  gradient of the shape function w.r.t. parent cs
%    COORD - nodal coordinate matrix (row oriented)
%    SDIM - spacial dimension (default is # cols in coord)
%    EDIM - dimension of parent cs (default is #cols dNdxi)
%
% function [J,JMAT]=ELEMENT_JACOBIAN(ETYPE,XI,COORD,SDIM)
%
% Perfoms the same function but takes the element type and a coordinate in
% the parent elemetn space to compute the Jacobian
%
%    ETYPE - the element type sstring ie 'Tria3'
%    DNDXI -  gradient of the shape function w.r.t. parent cs
%    COORD - nodal coordinate matrix (row oriented)
%    SDIM - spacial dimension (default is # cols in coord)
% 
% This is part ot FEMLAB
% Jack Chessa, jfchessa@utep.edu
%

if ( isstr(dNdxi) )
    
    etype=dNdxi;
    xi=coord;
    coord=sdim;
    dNdxi=dshape(etype,coord);
    
    if (nargin==4)
        sdim=edim;
    end
    
    edim=size(dNdxi,2);
    
else
    if ( nargin<3 )
        sdim=size(coord,2);
    end
    
    if (nargin<4)
        edim=size(dNdxi,2);
    end
end

jmat=zeros(sdim);
jmat(1:sdim,1:edim)=coord(:,1:sdim)'*dNdxi(:,1:edim);

if ( sdim~=edim )

    if ( sdim==2 & edim==1 )
        jmat(1,2) = -jmat(2,1);
		jmat(2,2) = jmat(1,1);
    elseif ( sdim==3 & edim==1 )
        jmat(:,2) = cross( [0;0;1], jmat(:,1) );
        if ( norm(jmat(:,2))==0 )
            jmat(:,2) = cross( [0;1;0], jmat(:,1) );
        end
        jmat(:,2) = jmat(:,2)/norm(jmat(:,2));
        jmat(:,3) = cross( jmat(:,1), jmat(:,2) );
        jmat(:,3) = jmat(:,3)/norm(jmat(:,3));
    elseif ( sdim==3 & edim==2 )
        jmat(:,3) = cross( jmat(:,1), jmat(:,2) );
        jmat(:,3) = jmat(:,3)/norm(jmat(:,3));
    end

end


jac=det(jmat);
% if ( edim~=sdim )
%     jac=jac^(edim/sdim);
% end
