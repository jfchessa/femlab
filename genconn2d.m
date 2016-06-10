function element=genconn2d(node_pattern,num_u,num_v,inc_u,inc_v)

% function ELEMENT=GENCONN2d(NODE_PATTERN,NUM_U,NUM_V,INC_U,INC_V)
%
% creates a connectivity list

if ( nargin < 3 )
   disp(['Not enough parameters specified for GENCONN2D function'])
end

if ( nargin < 5 )
    inc_v=2;
end
if ( nargin < 4 )
    inc_u=1;
end

inc=[zeros(1,size(node_pattern,2))];
e=1;
element=zeros(num_u*num_v,size(node_pattern,2));

for row=1:num_v
   for col=1:num_u
      element(e,:)=node_pattern+inc;
      inc=inc+inc_u;
      e=e+1;
   end
   inc=inc+inc_v-inc_u;
end
