% Program to read wind data from WRF netcdf file for assimilation with
% observations for 2014
%
% Read wind data;long, lat, wind speeds, wind stresses, and times; from 
% NetCDF file and save them in a matfile
%
% Then using Assimilate_WRF_Wind.m I calculate assimilation coafficients 
% and modified the winds 
%
% Next using "Calc_AF_along_sound.m" and "Assimilate_WRF_Wind_Field.m" 
% calculate wind field for LIS 
%
% finaly make a new wind forcing netcdf file using
% Read_WRF_make_new_netcdf_forcing.m and write_WRF_new_forcing.m
% 
% Oct 2109
% Amin Ilia - University of Connecticut, Groton, CT % 

%% Open and read WRF NetCDF file

clear all

addpath F:\Backup\WRF   % path for WRF forcing file

old_wrf='wrf_hnd_2014.nc';   % name of WRF forcing file

nc = netcdf.open(old_wrf, 'NOWRITE');

[numdims, numvars, numglobatts, unlimdimID] = netcdf.inq(nc);  % read the number of vars and dims in the netcdf file

for ii = 1:numvars  %% only read long, lat, wind speeds, wind stresses, and times required for assimilation

    % Find name of the current variable.
    [varname, xtype, varDimIDs, varAtts] = netcdf.inqVar(nc, ii - 1);  % read variables names
    
    % check the names for saving
    if (isequal(varname,'XLAT') || isequal(varname,'XLONG') ...        
         ||  isequal(varname,'U10') || isequal(varname,'V10') ...
         ||  isequal(varname,'Stress_U') || isequal(varname,'Stress_V') ...
         || isequal(varname,'Times')  )
   
        varid = netcdf.inqVarID(nc, varname);
    
        variable=netcdf.getVar(nc,varid);  % get data
    
        eval(['WRF.',varname,'=variable;']); % save in structure
    
    end

end

%% save wind data in a mat file

save('WRF','WRF','-v7.3') % save the structure WRF in file