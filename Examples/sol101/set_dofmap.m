function [dofMap,ngdof] = set_dofmap(solidlElem,shellElem,barElem,rodElem)

dofmap=[];
for e=1:length(solidElem)
    dofmap( solidElem(e).conn, 1:3 ) = 1;
end

for e=1:length(shellElem)
    dofmap( shellElem(e).conn, 1:6 ) = 1;
end

for e=1:length(barElem)
    dofmap( barElem(e).conn, 1:6 ) = 1;
end

for e=1:length(rodElem)
    dofmap( rodElem(e).conn, 1:3 ) = 1;
end

adofs=find(dofmap);
ngdof=length(adof);
dofmap(adofs)=1:ngdof;

end

