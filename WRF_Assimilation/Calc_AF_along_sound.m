% with a uniform function from -73.8 until -73.1 then sudenly decrease
% until -72.7 with quaric function and lat 40.8 to 41.3
%% Program to calculate assimilation factors from west to CLIS
% Oct, 2019
% Amin Ilia

%%
clear all

addpath C:\Users\amin\Desktop\Thesis\FVCOM\Data\CLIS

load('assi_coaf_v3.mat')


lonm=-73.8:0.0125:-72.7;
cfxl=ones(length(cfx),length(lonm));
cfyl=ones(length(cfx),length(lonm));
for i=1:57
    cfxl(:,i)=cfx;
    cfyl(:,i)=cfy;
end

ix=find(cfx~=1);
iy=find(cfy~=1);

for i=58:length(lonm)
     cfxl(:,i)=(2.8*(lonm(i)-lonm(57))-1).^4*(cfx-1)+1;
     cfyl(:,i)=(2.8*(lonm(i)-lonm(57))-1).^4*(cfy-1)+1;
end

save('cf_lon.mat','cfxl','cfyl','mtime','lonm')


