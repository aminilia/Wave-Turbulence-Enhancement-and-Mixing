%% The program for calculating wind modification coefficients for western LIS
%
% WRF underestimates wind speed at western LIS by 75%
% Here calculate modifications coefficients for wind at WLIS buoy
% The eastern, western, northern, and southern winds assimilated with
% different coefficient as the underestimation coefficient of WRF are
% directionally dependent. Only wind speeds larger than 3m/s are modified
%
% Oct 2019
% Amin Ilia - University of Connecticut, Groton, CT% 
% 
clear
close all

%% addpath and load data

addpath C:\Users\amin\Desktop\Thesis\FVCOM\Model\2014\Output_2014
addpath C:\Users\amin\Desktop\Thesis\Glider_Data\Air-Sea\air_sea
addpath C:\Users\amin\Desktop\Thesis\Data\COARE\matlab3_0\vectorized

load('model_2014_wholeyear.mat')
load('model_2014_temp_sal_wholeyear.mat')
load('model_2014_uvw_wholeyear.mat')
load('model_2014_stocks_wholeyear.mat')

%% Find Assimilate coefficients for whole year
addpath C:\Users\amin\Desktop\Thesis\FVCOM\Data\WLIS\Wind
load('wlis_wx_2006-2017.mat')

addpath C:\Users\amin\Desktop\Thesis\FVCOM\Data\CLIS
load('clis_wx.mat')

ws_w=convwind(wlis2014.windSpdKts);
ws_c=convwind(clis2014.windSpdKts);

ws_w=movmean(ws_w,4);   % One hour average
ws_c=movmean(ws_c,4);   % One hour average

wx_w=-ws_w.*sind(wlis2014.windDirM);
wy_w=-ws_w.*cosd(wlis2014.windDirM);

wx_c=-ws_c.*sind(clis2014.windDirM);
wy_c=-ws_c.*cosd(clis2014.windDirM);

wsq=sqrt(wxq.^2+wyq.^2);
wdeg=atan2d(wyq,wxq);

% WLIS
iw_m=find(wsq(:,6)>3.0);
twb=4.5;
iw_b=find(ws_w>twb);
c_w=mean(ws_w(iw_b))/mean(wsq(iw_m,6));

iw_me=find(wsq(:,6)>3.0 & wxq(:,6)>=0);
iw_mw=find(wsq(:,6)>3.0 & wxq(:,6)<0);
iw_mn=find(wsq(:,6)>3.0 & wyq(:,6)>=0);
iw_ms=find(wsq(:,6)>3.0 & wyq(:,6)<0);

iw_be=find(ws_w>twb & wx_w>=0);
iw_bw=find(ws_w>twb & wx_w<0);
iw_bn=find(ws_w>twb & wy_w>=0);
iw_bs=find(ws_w>twb & wy_w<0);

ce_w=mean(abs(wx_w(iw_be))/mean(abs(wxq(iw_me,6))));
cw_w=mean(abs(wy_w(iw_bw))/mean(abs(wyq(iw_mw,6))));
cn_w=mean(abs(wx_w(iw_bn))/mean(abs(wxq(iw_mn,6))));
cs_w=mean((wy_w(iw_bs))/mean((wyq(iw_ms,6))));

cx_w=mean(abs(wx_w(iw_b))/mean(abs(wxq(iw_m,6))));
cy_w=mean(abs(wy_w(iw_b))/mean(abs(wyq(iw_m,6))));


%% Correcting wind components
% for the times that WLIS wind data available, the coefficient is directly
% calculated from observation. The events that observation is not available
% the coefficients are estimated from yearly coefficients. 

wx_mod=wxq(:,6);    % assimilated winds x
wy_mod=wyq(:,6);    % assimilated winds y

cfx=ones(size(wx_mod),'like', wx_mod);  % assimilation factor array x
cfy=ones(size(wy_mod),'like', wy_mod);  % assimilation factor array y

cal=1;  % calibratcfyion factor when there is no observation
trmax=2.0; % treshold for max assimalation factor

for i=1:length(iw_m)
    
    dist=abs((wlis2014.EST+5/24)-mtime(iw_m(i)));   
    ii=find(dist==min(dist));
    
    if min(dist)<0.5
        if ((wx_w(ii)/wx_mod(iw_m(i)))<0)
            if (abs(wx_mod(iw_m(i)))>2)
                continue
            else      % if the wind on x is less than 2 the only calculate cfy
            	cfy(iw_m(i))=wy_w(ii)/wy_mod(iw_m(i));
                if (wy_mod(iw_m(i))>0 && cfy(iw_m(i))>trmax*cn_w)
                    cfy(iw_m(i))=trmax*cn_w;
                end
                if (wy_mod(iw_m(i))<0 && cfy(iw_m(i))>trmax*cs_w)
                    cfy(iw_m(i))=trmax*cs_w;
                end
                wy_mod(iw_m(i))=wy_mod(iw_m(i))*cfy(iw_m(i));
                continue
            end
        end
        if ((wy_w(ii)/wy_mod(iw_m(i)))<0)
            if (abs(wy_mod(iw_m(i)))>2)
            	continue
            else      % if the wind on y is less than 2 the only calculate cfx
            	cfx(iw_m(i))=wx_w(ii)/wx_mod(iw_m(i));
                if (wx_mod(iw_m(i))>0 && cfx(iw_m(i))>trmax*ce_w)
                	cfx(iw_m(i))=trmax*ce_w;
                end
                if (wx_mod(iw_m(i))<0 && cfx(iw_m(i))>trmax*cw_w)
                	cfx(iw_m(i))=trmax*cw_w;
                end
                wx_mod(iw_m(i))=wx_mod(iw_m(i))*cfx(iw_m(i));
                continue
            end
        end
        
        
        cfx(iw_m(i))=wx_w(ii)/wx_mod(iw_m(i));
        if (wx_mod(iw_m(i))>0 && cfx(iw_m(i))>trmax*ce_w)
            cfx(iw_m(i))=trmax*ce_w;
        end
        if (wx_mod(iw_m(i))<0 && cfx(iw_m(i))>trmax*cw_w)
             cfx(iw_m(i))=trmax*cw_w;
        end
        wx_mod(iw_m(i))=wx_mod(iw_m(i))*cfx(iw_m(i));
        
        cfy(iw_m(i))=wy_w(ii)/wy_mod(iw_m(i));
        if (wy_mod(iw_m(i))>0 && cfy(iw_m(i))>trmax*cn_w)
            cfy(iw_m(i))=trmax*cn_w;
        end
        if (wy_mod(iw_m(i))<0 && cfy(iw_m(i))>trmax*cs_w)
             cfy(iw_m(i))=trmax*cs_w;
        end
        wy_mod(iw_m(i))=wy_mod(iw_m(i))*cfy(iw_m(i));
        
    else
        if(wxq(iw_m(i),6)>=0)
            wx_mod(iw_m(i))=wx_mod(iw_m(i))*ce_w*cal;
            cfx(iw_m(i))=ce_w*cal;
        else
            wx_mod(iw_m(i))=wx_mod(iw_m(i))*cw_w*cal;
            cfx(iw_m(i))=cw_w*cal;
        end
        if(wyq(iw_m(i),6)>=0)
            wy_mod(iw_m(i))=wy_mod(iw_m(i))*cn_w*cal;
            cfy(iw_m(i))=cn_w*cal;
        else
            wy_mod(iw_m(i))=wy_mod(iw_m(i))*cs_w*cal;
            cfy(iw_m(i))=cs_w*cal;
        end
        
    end
    
end

cfx(cfx<1)=1;
cfy(cfy<1)=1;

save('assi_coaf_v3.mat','cfx','cfy','wx_mod','wy_mod','mtime')


%% plots to compare before and after assimilation

s_time=datenum(2014,1,1);
e_time=datenum(2015,1,1);

figure
subplot(2,1,1)
plot(wlis2014.EST+5/24,ws_w); hold on
plot(mtime,sqrt(wx_mod.^2+wy_mod.^2))
xlim([s_time e_time]);
ylabel('Wind Speed (m/s)')
grid on;hold on
datetick('x','mmm dd','keeplimits','keepticks')
title('WLIS')

subplot(2,1,2)
plot(clis2014.EST+5/24,ws_c);  hold on
plot(mtime,sqrt(wxq(:,7).^2+wyq(:,7).^2));
xlim([s_time e_time]);
ylabel('Wind Speed (m/s)')
grid on;hold on
datetick('x','mmm dd','keeplimits','keepticks')
title('CLIS')
print('-dpng','-r400','LIS_WLIS_CLIS_2014_wind_speeds_mvm_mod1.png')


figure
subplot(2,1,1)
plot(wlis2014.EST+5/24,ws_w); hold on
plot(mtime,sqrt(wxq(:,6).^2+wyq(:,6).^2));
xlim([s_time e_time]);
ylabel('Wind Speed (m/s)')
grid on;hold on
datetick('x','mmm dd','keeplimits','keepticks')
title('WLIS wind speed after assimilation')

subplot(2,1,2)
plot(wlis2014.EST+5/24,ws_w);  hold on
plot(mtime,sqrt(wx_mod.^2+wy_mod.^2))
xlim([s_time e_time]);
ylabel('Wind Speed (m/s)')
grid on;hold on
datetick('x','mmm dd','keeplimits','keepticks')
title('WLIS wind speed after assimilation')
print('-dpng','-r400','LIS_WLIS_2014_wind_speeds_assimilation_AI_v3_nlo_t.png')


figure
subplot(2,1,1)
plot(wlis2014.EST+5/24,-ws_w.*sind(wlis2014.windDirM)); hold on
plot(mtime,wx_mod); hold on
plot(mtime,wxq(:,6)); hold on
ylabel('Wind u (m/s)')
xlim([s_time e_time]);
ylim([-15 15])
grid on;hold on
datetick('x','mmm dd','keeplimits','keepticks')
title('WLIS')

subplot(2,1,2)
plot(wlis2014.EST+5/24,-ws_w.*cosd(wlis2014.windDirM)); hold on
plot(mtime,wy_mod); hold on
plot(mtime,wyq(:,6))
ylabel('Wind v (m/s)')
xlim([s_time e_time]);
ylim([-15 15])
grid on;hold on
datetick('x','mmm dd','keeplimits','keepticks')
print('-dpng','-r400','LIS_WLIS_2014_wind_mvm_mod2_t.png')


figure
subplot(2,1,1)
plot(wlis2014.EST+5/24,-ws_w.*sind(wlis2014.windDirM)); hold on
plot(mtime,wxq(:,6)); hold on
ylabel('Wind u (m/s)')
xlim([s_time e_time]);
ylim([-15 15])
grid on;hold on
datetick('x','mmm dd','keeplimits','keepticks')
title('WLIS u before assimilation')

subplot(2,1,2)
plot(wlis2014.EST+5/24,-ws_w.*sind(wlis2014.windDirM)); hold on
plot(mtime,wx_mod); hold on
ylabel('Wind u (m/s)')
xlim([s_time e_time]);
ylim([-15 15])
grid on;hold on
datetick('x','mmm dd','keeplimits','keepticks')
title('WLIS u after assimilation')
print('-dpng','-r400','LIS_WLIS_2014_wind_u_assimilation_AI_v3_nlo_t.png')


figure
subplot(2,1,1)
plot(wlis2014.EST+5/24,-ws_w.*cosd(wlis2014.windDirM)); hold on
% plot(mtime,wy_mod); hold on
plot(mtime,wyq(:,6))
ylabel('Wind v (m/s)')
xlim([s_time e_time]);
ylim([-15 15])
grid on;hold on
datetick('x','mmm dd','keeplimits','keepticks')
title('WLIS v before assimilation')

subplot(2,1,2)
plot(wlis2014.EST+5/24,-ws_w.*cosd(wlis2014.windDirM)); hold on
plot(mtime,wy_mod); hold on
% plot(mtime,wyq(:,6))
ylabel('Wind v (m/s)')
xlim([s_time e_time]);
ylim([-15 15])
grid on;hold on
datetick('x','mmm dd','keeplimits','keepticks')
title('WLIS v after assimilation')
print('-dpng','-r400','LIS_WLIS_2014_wind_v_assimilation_AI_v3_nlo_t.png')


figure
subplot(2,1,1)
plot(clis2014.EST+5/24,-ws_c.*sind(clis2014.windDirM)); hold on
plot(mtime,wxq(:,7)); 
xlim([s_time e_time])
ylabel('Wind u (m/s)')
ylim([-15 15])
grid on;hold on
datetick('x','mmm dd','keeplimits','keepticks')
title('CLIS')

subplot(2,1,2)
plot(clis2014.EST+5/24,-ws_c.*cosd(clis2014.windDirM)); hold on
plot(mtime,wyq(:,7));
xlim([s_time e_time])
ylabel('Wind v (m/s)')
ylim([-15 15])
grid on;hold on
datetick('x','mmm dd','keeplimits','keepticks')
print('-dpng','-r400','LIS_CLIS_2014_wind_mvm_1hour_mod1.png')


%% functions

function ws=convwind(ws)  
% this function convert winds from knot at 3.5m to m/s at 10m
    ws=ws/1.946;                             %Wind speed (m/s)
    z0=0.0001;                               %Roughness length
    ws=ws*log(10/z0)/log(3.5/z0);            %Correcting Wind speed for 10m heigh based on log wind profile (law of the wall)
end

