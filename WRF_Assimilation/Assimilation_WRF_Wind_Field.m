%% Program to assimilate WRF wind field
% Only assimilate data from lon -73.8 to -72.7  and lat 40.8 to 41.3 in LIS
% Oct 2019, 
% Amin Ilia - University of Connecticut, Groton, CT% 


%% Load data
clear all

% load wind data
load('WRF.mat')

% load assimilation coafficient along LIS
load('cf_lon.mat')

%% Found points that assimilate

ap=(WRF.XLONG>=-73.8 & WRF.XLONG<=-72.655 & WRF.XLAT>=40.8 & WRF.XLAT<=41.3);


%% Assimilate Western Points
% Note: time steps are the same in WRF and FVCOM therefore no need for
% temporal interpolation.
% However the WRF extends 10 days into 2015 therefore the number of
% timesteps are different

WRFm=WRF;

[ia,ja]=find(ap==1);
for i=1:length(ia)
    [~,ip(i)]=min(abs(WRF.XLONG(ia(i),ja(i))-lonm)); % find the closest point in calculated assimlation coafficients
    
    for k=1:size(cfxl,1)   % time
        WRFm.U10(ia(i),ja(i),k)=WRFm.U10(ia(i),ja(i),k).*cfxl(k,ip(i));      % modify wind speed vector u
        WRFm.V10(ia(i),ja(i),k)=WRFm.V10(ia(i),ja(i),k).*cfyl(k,ip(i));      % modify wind speed vector v
        
        WRFm.Stress_U(ia(i),ja(i),k)=WRFm.Stress_U(ia(i),ja(i),k).*cfxl(k,ip(i)).^2;     % modify wind stress vector u
        WRFm.Stress_V(ia(i),ja(i),k)=WRFm.Stress_V(ia(i),ja(i),k).*cfyl(k,ip(i)).^2;     % modify wind stress vector v
    end
end

%% checking a couple of storm events visualy
ttc=datenum('08-Apr-2014 03:00:00');
[~,itt]=min(abs(ttc-mtime));

figure
subplot(2,1,1)
pcolor(WRFm.XLONG,WRFm.XLAT,WRFm.U10(:,:,itt)); shading flat
colorbar
xlim([-74 -71])
ylim([40 42])
subplot(2,1,2)
pcolor(WRF.XLONG,WRF.XLAT,WRF.U10(:,:,itt)); shading flat
colorbar
xlim([-74 -71])
ylim([40 42])


ttc=datenum('03-Nov-2014 03:00:00');
[~,itt]=min(abs(ttc-mtime));

figure
subplot(2,1,1)
pcolor(WRFm.XLONG,WRFm.XLAT,WRFm.U10(:,:,itt)); shading flat
colorbar
xlim([-74 -71])
ylim([40 42])
subplot(2,1,2)
pcolor(WRF.XLONG,WRF.XLAT,WRF.U10(:,:,itt)); shading flat
colorbar
xlim([-74 -71])
ylim([40 42])


%% save modified WRF

% make julian time for netcdf 
st=juliandate(datetime('2014-01-01 00:00:00'),'modifiedjuliandate');
et=juliandate(datetime('2015-01-08 01:00:00'),'modifiedjuliandate');
WRFm.time=st:1/24:et;
WRFm.time=WRFm.time';

save('WRFm.mat', 'WRFm' ,'-v7.3')