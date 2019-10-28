function write_WRF_new_forcing(WRF, filename)
% Write modified wind data out to WRF format netCDF forcing file.
% This program is modified version of Pierre Cazenave (Plymouth Marine
% Laboratory) to replace assimilated wind data with orginal data in nc file
%
% write_WRF_new_forcing(WRF, filename)
%
% DESCRIPTION:
%   Takes the given regularly gridded forcing data and writes out to a WRF
%   format netCDF file.
%
% INPUT:
%   WRF - struct with the following fields:
%       XLONG         : longitude, rectangular array (see MESHGRID).
%       XLAT          : latitude, rectangular array (see MESHGRID).
%       time          : Modified Julian Day times.
%       SLP           : surface pressure [mb]
%       Shortwave     : shortwave radiation (upward = negative) [W/m^{2}]
%       Longwave      : longwave radiation (upward = negative) [W/m^{2}]
%       Sensible      : shortwave radiation (upward = negative) [W/m^{2}]
%       Latent        : longwave radiation (upward = negative) [W/m^{2}]
%       Net_Heat      : net surface heat flux (ocean losing = negative) [W/m^{2}]
%       U10           : eastward wind velocity [m/s]
%       V10           : northward wind velocity [m/s]
%       Stress_U      : eastward wind stress at sea surface [Pa]
%       Stress_V      : northward wind stress at sea surface [Pa]
%       precipitation : precipitation (ocean losing = negative) [m/s]
%       Evaporation   : evaporation (ocean losing = negative) [m/s]
%       RH            : relative humidity [%]
%       T2            : temperature at 2m [K]
%   filename - Output netCDF file name.
%
% OUTPUT:
%   WRF format heating netCDF file.
%
% EXAMPLE USAGE:
%   new_wrf_file = '/path/to/output/casename_wnd_new.nc';
%   write_WRF_new_forcing(WRF, new_wrf_file);
%
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
%   Amin Ilia (University of Connecticut)
%
% Revision history:
%   2019-10-21 - the program "write_WRF_new_forcing.m" produced by Pierre
%   Cazenave (Plymouth Marine Laboratory) was modified to replace 
%   assimilated wind speed and stress vectors.
%   % Amin Ilia - University of Connecticut, Groton, CT %
%
%==========================================================================

assert(nargin == 2, 'Incorrect number of arguments')

subname = 'write_WRF_new_forcing';

global ftbverbose
if ftbverbose
    fprintf('\nbegin : %s \n', subname)
end

ntimes = numel(WRF.time);
[sgYr, sgMon, sgDay, sgHr, sgMin, sgSec] = mjulian2greg(WRF.time(1));
[egYr, egMon, egDay, egHr, egMin, egSec] = mjulian2greg(WRF.time(end));

x=WRF.XLONG;
y=WRF.XLAT;
% Make the range of lon -180 - 180.
if max(x) > 180
    x(x > 180) = x(x > 180) - 360;
end
% I've yet to come across a latitude range that isn't -90 - 90, but just in
% case.
if max(y) > 90
    y(y > 90) = y(y > 90) - 180;
end
[nwest_east, nsouth_north] = size(x);

%--------------------------------------------------------------------------
% Create the netCDF header for the FVCOM forcing file
%--------------------------------------------------------------------------

nc = netcdf.create(filename, 'clobber');

netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'type', 'FVCOM METEO FORCING FILE')
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'title', [filename, ' forcing'])
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'gauge', 'Met Office Unified Model forcing')
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'history', sprintf('File with %s from modified MATLAB fvcom-toolbox', subname))
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'source', 'wrf grid (structured) surface forcing')
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'START_DATE', datestr(datenum(sgYr, sgMon, sgDay, sgHr, sgMin, sgSec), 'yyyy-mm-dd HH:MM:SS'))
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'END_DATE', datestr(datenum(egYr, egMon, egDay, egHr, egMin, egSec), 'yyyy-mm-dd HH:MM:SS'))
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'file', filename)
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'institution', 'University of Connecticut, Amin Ilia')
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'Conventions', 'CF-1.0')

% Dimensions
sn_dimid=netcdf.defDim(nc, 'south_north', nsouth_north);
we_dimid=netcdf.defDim(nc, 'west_east', nwest_east);
datestrlen_dimid=netcdf.defDim(nc, 'DateStrLen', 19);
time_dimid=netcdf.defDim(nc, 'time', netcdf.getConstant('NC_UNLIMITED'));

% Space variables
y_varid = netcdf.defVar(nc, 'XLAT', 'NC_FLOAT', [we_dimid, sn_dimid]);
netcdf.putAtt(nc, y_varid, 'long_name', 'latitude');
netcdf.putAtt(nc, y_varid, 'description', 'LATITUDE, SOUTH IS NEGATIVE');
netcdf.putAtt(nc, y_varid, 'units', 'degrees_north');
netcdf.putAtt(nc, y_varid, 'type', 'data');

x_varid=netcdf.defVar(nc, 'XLONG', 'NC_FLOAT', [we_dimid, sn_dimid]);
netcdf.putAtt(nc, x_varid, 'long_name', 'longitude');
netcdf.putAtt(nc, x_varid, 'description', 'LONGITUDE, WEST IS NEGATIVE');
netcdf.putAtt(nc, x_varid, 'units', 'degrees_east');
netcdf.putAtt(nc, x_varid, 'type', 'data');

% Time variables
times_varid=netcdf.defVar(nc, 'Times', 'NC_CHAR', [datestrlen_dimid, time_dimid]);
netcdf.putAtt(nc, times_varid, 'description', 'GMT time');
% netcdf.putAtt(nc, times_varid, 'format', 'String: Calendar Time');
% netcdf.putAtt(nc, times_varid, 'time_zone', 'UTC');

nswrs_varid = netcdf.defVar(nc, 'Shortwave', 'NC_FLOAT', [we_dimid, sn_dimid, time_dimid]);
netcdf.putAtt(nc, nswrs_varid, 'description', 'Shortwave, upward is negative');
netcdf.putAtt(nc, nswrs_varid, 'units', 'W m-2');
netcdf.putAtt(nc, nswrs_varid, 'grid', 'wrf_grid');
netcdf.putAtt(nc, nswrs_varid, 'coordinates', 'XLONG XLAT');
netcdf.putAtt(nc, nswrs_varid, 'type', 'data');

nlwrs_varid = netcdf.defVar(nc, 'Longwave', 'NC_FLOAT', [we_dimid, sn_dimid, time_dimid]);
netcdf.putAtt(nc, nlwrs_varid, 'description', 'Longwave, upward is negative');
netcdf.putAtt(nc, nlwrs_varid, 'units', 'W m-2');
netcdf.putAtt(nc, nlwrs_varid, 'grid', 'wrf_grid');
netcdf.putAtt(nc, nlwrs_varid, 'coordinates', 'XLONG XLAT');
netcdf.putAtt(nc, nlwrs_varid, 'type', 'data');

sensible_varid = netcdf.defVar(nc, 'Sensible', 'NC_FLOAT', [we_dimid, sn_dimid, time_dimid]);
netcdf.putAtt(nc, sensible_varid, 'description', 'Sensible Heat flux, upward is negative');
netcdf.putAtt(nc, sensible_varid, 'units', 'W m-2');
netcdf.putAtt(nc, sensible_varid, 'grid', 'wrf_grid');
netcdf.putAtt(nc, sensible_varid, 'coordinates', 'XLONG XLAT');
netcdf.putAtt(nc, sensible_varid, 'type', 'data');

latent_varid = netcdf.defVar(nc, 'Latent', 'NC_FLOAT', [we_dimid, sn_dimid, time_dimid]);
netcdf.putAtt(nc, latent_varid, 'description', 'Latent Heat flux, upward is negative');
netcdf.putAtt(nc, latent_varid, 'units', 'W m-2');
netcdf.putAtt(nc, latent_varid, 'grid', 'wrf_grid');
netcdf.putAtt(nc, latent_varid, 'coordinates', 'XLONG XLAT');
netcdf.putAtt(nc, latent_varid, 'type', 'data');

nshf_varid = netcdf.defVar(nc, 'Net_Heat', 'NC_FLOAT', [we_dimid, sn_dimid, time_dimid]);
netcdf.putAtt(nc, nshf_varid, 'description', 'Sum of shortwave, longwave, sensible and latent heat fluxes, ocean lose heat is negative');
netcdf.putAtt(nc, nshf_varid, 'units', 'W m-2');
netcdf.putAtt(nc, nshf_varid, 'grid', 'wrf_grid');
netcdf.putAtt(nc, nshf_varid, 'coordinates', 'XLONG XLAT');
netcdf.putAtt(nc, nshf_varid, 'type', 'data');

u10_varid = netcdf.defVar(nc, 'U10', 'NC_FLOAT', [we_dimid, sn_dimid, time_dimid]);
netcdf.putAtt(nc, u10_varid, 'long_name', 'Eastward Wind Velocity');
netcdf.putAtt(nc, u10_varid, 'description', 'U at 10 M');
netcdf.putAtt(nc, u10_varid, 'units', 'm s-1');
netcdf.putAtt(nc, u10_varid, 'grid', 'wrf_grid');
netcdf.putAtt(nc, u10_varid, 'coordinates', 'XLONG XLAT');
netcdf.putAtt(nc, u10_varid, 'type', 'data');

v10_varid = netcdf.defVar(nc, 'V10', 'NC_FLOAT', [we_dimid, sn_dimid, time_dimid]);
netcdf.putAtt(nc, v10_varid, 'long_name', 'Northward Wind Velocity');
netcdf.putAtt(nc, v10_varid, 'description', 'V at 10 M');
netcdf.putAtt(nc, v10_varid, 'units', 'm s-1');
netcdf.putAtt(nc, v10_varid, 'grid', 'wrf_grid');
netcdf.putAtt(nc, v10_varid, 'coordinates', 'XLONG XLAT');
netcdf.putAtt(nc, v10_varid, 'type', 'data');

Stress_U_varid = netcdf.defVar(nc, 'Stress_U', 'NC_FLOAT', [we_dimid, sn_dimid, time_dimid]);
netcdf.putAtt(nc, Stress_U_varid, 'long_name', 'Eastward Wind Stress');
netcdf.putAtt(nc, Stress_U_varid, 'description', 'U Wind stress at sea surface, westward is negative');
netcdf.putAtt(nc, Stress_U_varid, 'units', 'Pa');
netcdf.putAtt(nc, Stress_U_varid, 'grid', 'wrf_grid');
netcdf.putAtt(nc, Stress_U_varid, 'coordinates', 'XLONG XLAT');
netcdf.putAtt(nc, Stress_U_varid, 'type', 'data');

Stress_V_varid = netcdf.defVar(nc, 'Stress_V', 'NC_FLOAT', [we_dimid, sn_dimid, time_dimid]);
netcdf.putAtt(nc, Stress_V_varid, 'long_name', 'Northward Wind Stress');
netcdf.putAtt(nc, Stress_V_varid, 'description', 'V Wind stress at sea surface, westward is negative');
netcdf.putAtt(nc, Stress_V_varid, 'units', 'Pa');
netcdf.putAtt(nc, Stress_V_varid, 'grid', 'wrf_grid');
netcdf.putAtt(nc, Stress_V_varid, 'coordinates', 'XLONG XLAT');
netcdf.putAtt(nc, Stress_V_varid, 'type', 'data');

prate_varid = netcdf.defVar(nc, 'Precipitation', 'NC_FLOAT', [we_dimid, sn_dimid, time_dimid]);
netcdf.putAtt(nc, prate_varid, 'long_name', 'Precipitation');
netcdf.putAtt(nc, prate_varid, 'description', 'Precipitation, ocean lose water is negative');
netcdf.putAtt(nc, prate_varid, 'units', 'm s-1');
netcdf.putAtt(nc, prate_varid, 'grid', 'wrf_grid');
netcdf.putAtt(nc, prate_varid, 'coordinates', 'XLONG XLAT');
netcdf.putAtt(nc, prate_varid, 'type', 'data');

evap_varid = netcdf.defVar(nc, 'Evaporation', 'NC_FLOAT', [we_dimid, sn_dimid, time_dimid]);
netcdf.putAtt(nc, evap_varid, 'long_name', 'Evaporation');
netcdf.putAtt(nc, evap_varid, 'description', 'Evaporation, ocean lose water is negative');
netcdf.putAtt(nc, evap_varid, 'units', 'm s-1');
netcdf.putAtt(nc, evap_varid, 'grid', 'wrf_grid');
netcdf.putAtt(nc, evap_varid, 'coordinates', 'XLONG XLAT');
netcdf.putAtt(nc, evap_varid, 'type', 'data');

pres_varid = netcdf.defVar(nc, 'SLP', 'NC_FLOAT', [we_dimid, sn_dimid, time_dimid]);
netcdf.putAtt(nc, pres_varid, 'long_name', 'Air Pressure');
netcdf.putAtt(nc, pres_varid, 'description', 'Sea level pressure only for ocean');
netcdf.putAtt(nc, pres_varid, 'units', 'mb');
netcdf.putAtt(nc, pres_varid, 'grid', 'wrf_grid');
netcdf.putAtt(nc, pres_varid, 'coordinates', 'XLONG XLAT');
netcdf.putAtt(nc, pres_varid, 'type', 'data');

rh_varid = netcdf.defVar(nc, 'RH', 'NC_FLOAT', [we_dimid, sn_dimid, time_dimid]);
netcdf.putAtt(nc, rh_varid, 'long_name', 'Relative humidity');
netcdf.putAtt(nc, rh_varid, 'description', 'Relative humidity (generated from Qv)');
netcdf.putAtt(nc, rh_varid, 'units', '%');
netcdf.putAtt(nc, rh_varid, 'grid', 'wrf_grid');
netcdf.putAtt(nc, rh_varid, 'coordinates', 'XLONG XLAT');
netcdf.putAtt(nc, rh_varid, 'type', 'data');

air_varid = netcdf.defVar(nc, 'T2', 'NC_FLOAT', [we_dimid, sn_dimid, time_dimid]);
netcdf.putAtt(nc, air_varid, 'long_name', 'Air Temperature at 2m');
netcdf.putAtt(nc, air_varid, 'description', 'Bulk air temperature');
netcdf.putAtt(nc, air_varid, 'units', 'K');
netcdf.putAtt(nc, air_varid, 'grid', 'wrf_grid');
netcdf.putAtt(nc, air_varid, 'coordinates', 'XLONG XLAT');
netcdf.putAtt(nc, air_varid, 'type', 'data');


% End definitions
netcdf.endDef(nc);

%--------------------------------------------------------------------------
% Put the data in the netCDF file.
%--------------------------------------------------------------------------

% Build the Times string and output to netCDF.
nStringOut = char();
[nYr, nMon, nDay, nHour, nMin, nSec] = mjulian2greg(WRF.time);
for tt = 1:ntimes
    nDate = [nYr(tt), nMon(tt), nDay(tt), nHour(tt), nMin(tt), nSec(tt)];
    nStringOut = [nStringOut, sprintf('%04i/%02i/%02i %02i:%02i:%02i', nDate)];
end
netcdf.putVar(nc, times_varid, [0, 0], [19, ntimes], nStringOut);
% And the rest...
netcdf.putVar(nc, x_varid, x);
netcdf.putVar(nc, y_varid, y);   
netcdf.putVar(nc, air_varid, [0, 0, 0], [nwest_east, nsouth_north, ntimes], WRF.T2);
netcdf.putVar(nc, nswrs_varid, [0, 0, 0], [nwest_east, nsouth_north, ntimes], WRF.Shortwave);
netcdf.putVar(nc, nlwrs_varid, [0, 0, 0], [nwest_east, nsouth_north, ntimes], WRF.Longwave);
netcdf.putVar(nc, sensible_varid, [0, 0, 0], [nwest_east, nsouth_north, ntimes], WRF.Sensible);
netcdf.putVar(nc, latent_varid, [0, 0, 0], [nwest_east, nsouth_north, ntimes], WRF.Latent);
netcdf.putVar(nc, nshf_varid, [0, 0, 0], [nwest_east, nsouth_north, ntimes], WRF.Net_Heat);
netcdf.putVar(nc, u10_varid, [0, 0, 0], [nwest_east, nsouth_north, ntimes], WRF.U10);
netcdf.putVar(nc, v10_varid, [0, 0, 0], [nwest_east, nsouth_north, ntimes], WRF.V10);
netcdf.putVar(nc, Stress_U_varid, [0, 0, 0], [nwest_east, nsouth_north, ntimes], WRF.Stress_U);
netcdf.putVar(nc, Stress_V_varid, [0, 0, 0], [nwest_east, nsouth_north, ntimes], WRF.Stress_V);
netcdf.putVar(nc, prate_varid, [0, 0, 0], [nwest_east, nsouth_north, ntimes], WRF.SLP);
netcdf.putVar(nc, evap_varid, [0, 0, 0], [nwest_east, nsouth_north, ntimes], WRF.Evaporation);
netcdf.putVar(nc, pres_varid, [0, 0, 0], [nwest_east, nsouth_north, ntimes], WRF.Precipitation);
netcdf.putVar(nc, rh_varid, [0, 0, 0], [nwest_east, nsouth_north, ntimes], WRF.RH);



% Close the netCDF file
netcdf.close(nc);


fprintf('end   : %s \n', subname)
end