function element=genconn3d(node_pattern,num_u,num_v,num_w,inc_u,inc_v,inc_w)

% function ELEMENT=GENCONN3d(NODE_PATTERN,NUM_U,NUM_V,NUM_W,
%               INC_U,INC_V,INC_W)
%
% creates a connectivity list

if ( nargin < 4 )
   disp(['Not enough parameters specified for GENCONN3D function'])
end

if ( nargin < 7 )
    inc_w=num_u+3;
end
if ( nargin < 6 )
    inc_v=2;
end
if ( nargin < 5 )
    inc_u=1;
end


inc=[zeros(1,size(node_pattern,2))];
e=1;
element=zeros(num_u*num_v*num_w,size(node_pattern,2));

for k=1:num_w
    for j=1:num_v
        for i=1:num_u
            element(e,:)=node_pattern+inc;
            inc=inc+inc_u;
            e=e+1;
        end
        inc=inc+inc_v-inc_u;
    end
    inc=inc+inc_w-inc_v;
end
