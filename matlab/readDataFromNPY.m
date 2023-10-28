%% Load rawData3D
clear; close all;

addpath(genpath('.\npy-matlab'));
addpath(genpath('.\data'));
rawData = readNPY('B_008.npy');

frame_num   =       96;
chirp_num   =       64;
adcnum      =       32;
tx_num      =       2;
rx_num      =       4;
downsample  =       2;

adc_data = permute(rawData, [1, 2, 4, 3]); 
adc_data = downsample*adc_data(:,:,:,1:downsample:end);

shape   =   size(adc_data);    
w       =   reshape(hamming(shape(4)),1,1,1,[]);%generate window
data = adc_data.*w;
radar_cube = fft(data,64,4);
radar_mean = mean(radar_cube,2);
radar_cube = radar_cube-radar_mean;


fft1d_in = permute(radar_cube(:,:,[8,7,6,5,4,3,2,1],:),[1, 2, 4, 3]);

shape           =   size(fft1d_in);          
w               =   reshape(hamming(shape(4)),1,1,1,[]);%generate window
fft1d_in        =   fft1d_in.*w;
range_azimuth   =   fft(fft1d_in, 64, 4);
range_azimuth   =   fftshift(squeeze(sum(range_azimuth(:,1:3,:,:),2)),3);
range_azimuth   =   flip(range_azimuth,3);
range_azimuth   =   flip(range_azimuth,2);
range_azimuth   =   abs(range_azimuth);


set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultTextFontName', 'Times New Roman');

fig = figure(1);

set(fig, 'Color', 'w');
pos = fig.Position;

pos(1) = 200;
pos(2) = 200;
pos(3) = 320;
pos(4) = 320;


fig.Position = pos;

ax1 = axes('Parent', fig, 'Units', 'pixels', 'Position', [40,40, 240, 240]);

title(ax1, 'Range-Angle Feature', 'FontSize', 16)

h_img = image(ax1,NaN);
axis(ax1, 'equal');
xlim(ax1, [1 64])
ylim(ax1, [1 64])

xlabel(ax1, 'angle', 'FontSize', 14);
ylabel(ax1, 'range', 'FontSize', 14);

xticks(ax1, []);
yticks(ax1, []);


ax1.XAxis.TickLength = [0 0];
ax1.YAxis.TickLength = [0 0];

box(ax1, 'on');

for i=1:frame_num
    rai_abs = squeeze(range_azimuth(i,:,:));

    set(h_img, 'CData', rai_abs/max(max(rai_abs))*255);

    cmap = colormap(parula);
    pause(0.05)
end