% Program to read all data from orginal WRF, replace wind speeds and 
% stresses with assimilated ones for 2014.
% This progtram uses write_WRF_new_forcing.m which is a modified version of
% preprocessing toolbox in FVCOM.
% Oct 2109
% Amin Ilia - University of Connecticut, Groton, CT % 

%% 

clear all

addpath /d1/amin/LIS/2014/input   % path for WRF forcing file

old_wrf='wrf_hnd_2014.nc';   % name of WRF forcing file

out_wrf=[old_wrf(1:end-3),'_Assimilated.nc'];

%% Read data in orginal WRF file

nc = netcdf.open(old_wrf, 'NOWRITE');

[numdims, numvars, numglobatts, unlimdimID] = netcdf.inq(nc);  % read the number of vars and dims in the netcdf file

for ii = 1:numvars  % only save long, lat, wind speeds and stresses required for assimilation

    % Find name of the current variable.
    [varname, xtype, varDimIDs, varAtts] = netcdf.inqVar(nc, ii - 1);  % read variables names
    
    varid = netcdf.inqVarID(nc, varname);
    
    variable=netcdf.getVar(nc,varid);  % get data
    
    eval(['WRF.',varname,'=variable;']);

end

%% Replace assimilated wind data with orginal one

load('WRFm.mat') % load modified wind field

WRF.time=WRFm.time;
WRF.U10=WRFm.U10;
WRF.V10=WRFm.V10;
WRF.U10=WRFm.U10;
WRF.V10=WRFm.V10;
WRF.Stress_U=WRFm.Stress_U;
WRF.Stress_V=WRFm.Stress_V;

%% write new netcdf forcing file 

write_WRF_new_forcing(WRF, out_wrf)