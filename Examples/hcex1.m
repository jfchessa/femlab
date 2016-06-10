% hcex1.m
%
% Example input for the heat_conduction.m code
%
clear

gmshFile='MeshFiles/hcex1.msh';

matlRegion=[14,15];
matlAlpha(matlRegion)=[30,10];

fixedTemperatureRegion=[12 ];
fixedTemperature=[100];

heatFluxRegion=[];
heatFlux=[];

convectionRegion=[13];
covectionCoefficient=[10];
convectionTemperature=[20]

heat_conduction