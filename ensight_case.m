function stat=ensight_case(jobname,geomfile,times,...
  scalarvar,vectorvar,tensorvar,scalarcel,vectorcel,tensorcel)
% function STAT=ENSIGHT_CASE(JOBNAME,GEOMFILE,TIMES,
%         SCALARVAR,VECTORVAR,TENSORVAR, SCALARCELL,VECTORCELL,TENSORCELL )
%
% Write an Ensight .case file
%
% JOBNAME - name of the .case file.  Do  not add the .case extension
% GEOMFILE - name of the ensight geometry file
% TIMES - a vector of the timesteps 
% SCALARVALS - a cell structure of the names of the scalar variables
% VECTORVALS - a cell structure of the names of the vector variables
% TENSORVALS - a cell structure of the names of the tensor variables
% SCALARVALS - a cell structure of the names of the scalar per element
% VECTORVALS - a cell structure of the names of the vector per element
% TENSORVALS - a cell structure of the names of the tensor per element
%
% If the TIMES variable is not empty
% then it will assume that the variables file names are named accordingly
%
%    JOBNAME****.VARNAME
%
% else it will assume the following naming convention
%
%    JOBNAME.VARNAME
%
% Note: to generate file names in the JOBNAME****.VARNAME format use 
%       ['JOBNAME','num2str(n,'%04i'),'.VARNAME']
%
%

fid=fopen([jobname,'.case'],'w');

if ( fid<0 )
    stat=1;
    return
else
    stat=0;
end

if ( length(times)>0 )
    ts='****.';
else
    ts='.';
end

fprintf(fid,'FORMAT\ntype: ensight gold\n');

fprintf(fid,'\nGEOMETRY\n');
fprintf(fid,'model:   %s\n',geomfile);

fprintf(fid,'\nVARIABLE\n');
for i=1:length(scalarvar)
    fprintf(fid,'scalar per node: %s\n',[jobname,ts,scalarvar{i}]);
end

if ( nargin>4 )
    for i=1:length(vectorvar)
        fprintf(fid,'vector per node: %s\n',[jobname,ts,vectorvar{i}]);
    end
end

if ( nargin>5 )
    for i=1:length(tensorvar)
        fprintf(fid,'tensor per node: %s\n',[jobname,ts,tensorvar{i}]);
    end
end
if ( nargin>6 )
    for i=1:length(scalarcel)
        fprintf(fid,'scalar per element: %s\n',[jobname,ts,scalarcel{i}]);
    end
end
if ( nargin>7 )
    for i=1:length(vectorcel)
        fprintf(fid,'vector per element: %s\n',[jobname,ts,vectorcel{i}]);
    end
end

if ( nargin>8 )
    for i=1:length(tensorcel)
        fprintf(fid,'tensor per element: %s\n',[jobname,ts,tensorcel{i}]);
    end
end

if ( length(times)>0 )
    fprintf(fid,'\nTIME\n');
    fprintf(fid,'time set: 1\n');
    fprintf(fid,'number of steps: %5i\n',length(times));
    fprintf(fid,'filename start number: 0\n');
    fprintf(fid,'filename increment: 1\n');
    fprintf(fid,'time values:\n');
    fprintf(fid,'%12.6f\n',times);
end

fclose(fid); 