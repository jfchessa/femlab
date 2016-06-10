function phiBar = phi_mod(node,element,phi)

% function phiBar = phi_mod(node,element,phi)

numnode=size(node,1);
numelem=size(element,1);

psi=zeros(numnode,1);
%***********************************************************************
% COLLECT ELEMENT SET T AND NODE SET Ne
counter=1;
flag=0;
for e = 1 : numelem
    
    sctr = element(e,:);
    phis=phi(sctr);
    
    for judgeT = 1 : size(phis,1)-1
        if phis(judgeT)*phis(judgeT+1) < 0
            setofT(counter)=e;
            flag=1;
        end
    end
    if phis(judgeT+1)*phis(1) < 0
        setofT(counter)=e;
        flag=1;
    end
    if flag == 1
        counter=counter+1;
        flag=0;
    end
end

setofNe=unique(element(setofT,:));

%***********************************************************************
%ENRICHED FUNCTION VALUES IN THE NODE SET Ne
psi(setofNe)=abs(phi(setofNe));

%***********************************************************************
% COLLECT THE ADJACENT ELEMENT SET AND NODE SET Npsi
counter=1;

for e = 1 : numelem
    
    flag=ismember(element(e,:),setofNe);
    if any(flag)
        if ~ismember(e, setofT)
            setofNadj(counter)=e;
            counter=counter+1;
        end
    end
end

setofNpsi=unique(element(setofNadj,:));
tempNset=intersect(setofNpsi,setofNe);
setofNpsi=setxor(setofNpsi,tempNset);
%***********************************************************************
%ENRICHED FUNCTION VALUES IN THE NODE SET Npsi
psi(setofNpsi)=abs(phi(setofNpsi));
%***********************************************************************
%SORT NODE SET Npsi BASED ON THE ENRICHED INITIAL VALUES ON THE NODES
[sorted,I]=sort(psi(setofNpsi));
setofNpsi=setofNpsi(I);
%***********************************************************************
%GET THE ENRICHED FUNCTION VALUES FOR THE NODES IN THE NODE SET Npsi

for n = 1 : size(setofNpsi,1)
    
    % CREATE THE NODE SET PJ FOR EVERY NODE IN NODE SET Npsi
    [setofJsupport,noUse]=find(element==setofNpsi(n)); %setofJsupport IS THE ELEMENT SET CONTAINING NODE J
    setofSupport=unique(element(setofJsupport,:)); %setofSupport IS THE SET CONTAINING ALL THE NODES IN setofJsupport
    tempNset=union(setofNe,setofNpsi);
    tempNset=intersect(tempNset,setofSupport);
    tempIndex=find(psi(tempNset)<psi(setofNpsi(n)));
    
    if tempIndex==[]
        continue;
    end
    
    setofPJ=tempNset(tempIndex);
    
    L=zeros(size(setofPJ,1),1);
    for nn=1 : size(setofPJ,1)
        R=node(setofPJ(nn),:)-node(setofNpsi(n),:);
        L(nn)=sqrt(R*R');
    end
    Lk=sum(1./L.^2);
    psi(setofNpsi(n))=sum(1/Lk*1./L.^2.*psi(setofPJ));
end

phiBar=psi.*sign(phi);

restOfNodes = setdiff( 1:numnode, setofNpsi );
restOfNodes = setdiff( restOfNodes, setofNe );

phiBar(restOfNodes)=sign(phi(restOfNodes))*mean( abs(phiBar(setofNpsi)) );

% clf;
% plot_field([node phiBar],element,'T3',phiBar);
% hold on
% plot_mesh([node phiBar],element,'T3','k-');