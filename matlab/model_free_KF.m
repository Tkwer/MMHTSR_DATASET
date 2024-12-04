
%% Load rawData3D
clear; close all;%清屏
dataName    =       'res1'; % Change only this line
%% Define parameters, update based on the scenario
rawData     =       load(dataName);
rawData     =       rawData.gene_features;

addpath(genpath('.\npy-matlab'));
rawData = readNPY('B_006.npy');

frame_num   =       96;
chirp_num   =       64;
adcnum      =       32;
tx_num      =       2;
rx_num      =       4;
downsample  =       2;

adc_data = permute(rawData, [1, 2, 4, 3]); 
adc_data = downsample*adc_data(:,:,:,1:downsample:end);


PWR_threshold = 15.0;
NOISE_threshold = 30;
MAX_power = 1.6412e+5;

shape   =   size(adc_data);    
w       =   reshape(hamming(shape(4)),1,1,1,[]);%generate window
data = adc_data.*w;
radar_cube = fft(data,64,4);
radar_mean = mean(radar_cube,2);
% radar_mean = reshape(radar_mean,96,64,8,1);
radar_cube = radar_cube-radar_mean;
radar_cube(abs(radar_cube)<NOISE_threshold)=0;


% # 计算一帧的能量
frame_PWR = log(sum(abs(radar_cube),[2 3 4]));

scaling_factor = zeros(1,frame_num);
% Find indices of values greater than PWR_threshold
idx = frame_PWR > PWR_threshold;
% Calculate factor for values greater than PWR_threshold
scaling_factor(idx)  = 1;
% Set factor to 0 for values less than or equal to PWR_threshold
scaling_factor(~idx) = 0;

fft1d_in = permute(radar_cube(:,:,[8,7,6,5,4,3,2,1],:),[1, 2, 4, 3]);

shape           =   size(fft1d_in);          
w               =   reshape(hamming(shape(4)),1,1,1,[]);%generate window
fft1d_in        =   fft1d_in.*w;
range_azimuth   =   fft(fft1d_in, 64, 4);
range_azimuth   =   fftshift(squeeze(sum(range_azimuth(:,1:3,:,:),2)),3);
range_azimuth   =   flip(range_azimuth,3);
range_azimuth   =   flip(range_azimuth,2);
% # 抑制脉冲能量 归一化能量
% rai_abs = abs(range_azimuth(i,1,:,:)));
range_azimuth = abs(range_azimuth)./MAX_power;

kq = ones(64,64);
queue_center =zeros(5,2);
queue_center1 =zeros(0,2);


global P
global Q
global R

global dk
global xc

dk = 0;
xc = 0;
P = eye(2);
Q = [1e-2,0;0,1e-2];
R = [0.05,0;0,0.005];
queue_track = zeros(1,2);


set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultTextFontName', 'Times New Roman');

fig = figure(1);
% 设置整个窗口背景颜色为纯白
set(fig, 'Color', 'w');
pos = fig.Position;


% 设置窗口大小为 800x600 像素
pos(1) = 200;
pos(2) = 200;
pos(3) = 1320;
pos(4) = 320;

% 设置新的窗口位置和大小
fig.Position = pos;

% 创建子图并设置固定大小值

ax1 = axes('Parent', fig, 'Units', 'pixels', 'Position', [80, 180,320, 100]);
ax2 = axes('Parent', fig, 'Units', 'pixels', 'Position', [80, 40, 320, 100]);
ax3 = axes('Parent', fig, 'Units', 'pixels', 'Position', [440,40, 240, 240]);
ax4 = axes('Parent', fig, 'Units', 'pixels', 'Position', [720,40, 240, 240]);
ax5 = axes('Parent', fig, 'Units', 'pixels', 'Position', [1000,40, 240, 240]);

%% 绘制x的局部高斯过程回归
hold(ax1, 'on') % 将绘图保持在同一个图中
title(ax1, 'Angle-Time', 'FontSize', 14)

%画真值
h10 = scatter(ax1, NaN, NaN, ...                % 数据点的x和y坐标
    25, ...                      % 数据点大小
    [0.5, 0.2, 0.6], ...          % 数据点填充颜色（自定义颜色）
    'filled', ...                 % 填充数据点
    'MarkerEdgeColor', [0.5, 0.2, 0.6], ...   % 数据点边缘颜色（绿色）
    'MarkerFaceAlpha', 0.2, ...   % 设置数据点填充颜色透明度
    'MarkerEdgeAlpha', 1);        % 设置数据点边缘颜色透明度

% %画局部预测的值
% h11 = scatter(ax1, NaN, NaN, ...                % 数据点的x和y坐标
%     25, ...                      % 数据点大小
%     [109/255, 173/255, 209/255], ...          % 数据点填充颜色（自定义颜色）
%     'filled', ...                 % 填充数据点
%     'MarkerEdgeColor', [109/255, 173/255, 209/255], ...   % 数据点边缘颜色（绿色）
%     'MarkerFaceAlpha', 0.4, ...   % 设置数据点填充颜色透明度
%     'MarkerEdgeAlpha', 1);        % 设置数据点边缘颜色透明度

%画滤波的值
h12 = scatter(ax1, NaN, NaN, ...                % 数据点的x和y坐标
    25, ...                      % 数据点大小
    [183/255, 034/255, 048/255], ...          % 数据点填充颜色（自定义颜色）
    'filled', ...                 % 填充数据点
    'MarkerEdgeColor', [183/255, 034/255, 048/255], ...   % 数据点边缘颜色（绿色）
    'MarkerFaceAlpha', 0.4, ...   % 设置数据点填充颜色透明度
    'MarkerEdgeAlpha', 1);        % 设置数据点边缘颜色透明度

h13 = fill(ax1, [NaN,NaN,NaN], [NaN,NaN,NaN] ,[183/255, 034/255, 048/255],'FaceAlpha', 0.5, 'EdgeAlpha', 0.5, 'EdgeColor', [183/255, 034/255, 048/255]);%画不确定度区间
xlim(ax1, [0 96])
ylim(ax1, [0 64])
% 横纵坐标显示
x1_label = xlabel(ax1, 't', 'FontSize', 12);
ylabel(ax1, 'angle', 'FontSize', 12);

% % 获取当前位置
% pos = x1_label.Position;
% % 向上移动X轴标签（增加第2个坐标值）
% pos(2) = pos(2) + 10; % 在这里，0.1是您要向上移动的距离，可以根据需要调整
% % 更新X轴标签的位置
% x1_label.Position = pos;

% 关闭x轴和y轴的刻度值
xticks(ax1, []);
yticks(ax1, []);

% 关闭x轴和y轴的刻度线
ax1.XAxis.TickLength = [0 0];
ax1.YAxis.TickLength = [0 0];

% 保留封闭的box
box(ax1, 'on');

% 添加图标
legend(ax1, 'Raw','Filtered')

hold(ax1, 'off');

%% 绘制y的局部高斯过程回归
hold(ax2, 'on') % 将绘图保持在同一个图中
title(ax2, 'Range-Time', 'FontSize', 14)
%画真值
h20 = scatter(ax2, NaN, NaN, ...                % 数据点的x和y坐标
    25, ...                      % 数据点大小
    [0.5, 0.2, 0.6], ...          % 数据点填充颜色（自定义颜色）
    'filled', ...                 % 填充数据点
    'MarkerEdgeColor', [0.5, 0.2, 0.6], ...   % 数据点边缘颜色（绿色）
    'MarkerFaceAlpha', 0.2, ...   % 设置数据点填充颜色透明度
    'MarkerEdgeAlpha', 1);        % 设置数据点边缘颜色透明度

% %画局部预测的值
% h21 = scatter(ax2, NaN, NaN, ...                % 数据点的x和y坐标
%     25, ...                      % 数据点大小
%     [109/255, 173/255, 209/255], ...          % 数据点填充颜色（自定义颜色）
%     'filled', ...                 % 填充数据点
%     'MarkerEdgeColor', [109/255, 173/255, 209/255], ...   % 数据点边缘颜色（绿色）
%     'MarkerFaceAlpha', 0.4, ...   % 设置数据点填充颜色透明度
%     'MarkerEdgeAlpha', 1);        % 设置数据点边缘颜色透明度

%画滤波后的值
h22 = scatter(ax2, NaN, NaN, ...                % 数据点的x和y坐标
    25, ...                      % 数据点大小
    [183/255, 034/255, 048/255], ...          % 数据点填充颜色（自定义颜色）
    'filled', ...                 % 填充数据点
    'MarkerEdgeColor', [183/255, 034/255, 048/255], ...   % 数据点边缘颜色（绿色）
    'MarkerFaceAlpha', 0.4, ...   % 设置数据点填充颜色透明度
    'MarkerEdgeAlpha', 1);        % 设置数据点边缘颜色透明度

h23 = fill(ax2, [NaN,NaN,NaN], [NaN,NaN,NaN] ,[183/255, 034/255, 048/255],'FaceAlpha', 0.5, 'EdgeAlpha', 0.5, 'EdgeColor', [183/255, 034/255, 048/255]);%画不确定度区间

xlim(ax2, [0 96])
ylim(ax2, [0 64])

% 横纵坐标显示
x_label = xlabel(ax2, 't', 'FontSize', 12);
ylabel(ax2, 'range', 'FontSize', 12);

% % 获取当前位置
% pos = x_label.Position;
% % 向上移动X轴标签（增加第2个坐标值）
% pos(2) = pos(2) + 0.2; % 在这里，0.1是您要向上移动的距离，可以根据需要调整
% % 更新X轴标签的位置
% x_label.Position = pos;

% 关闭x轴和y轴的刻度值
xticks(ax2, []);
yticks(ax2, []);

% 关闭x轴和y轴的刻度线
ax2.XAxis.TickLength = [0 0];
ax2.YAxis.TickLength = [0 0];

% 保留封闭的box
box(ax2, 'on');

% 添加图标
legend(ax2, 'Raw','Filtered')
hold(ax2, 'off');

%% 绘制xy联合局部高斯过程回归
hold(ax3, 'on') % 将绘图保持在同一个图中
title(ax3, 'Angle-Range', 'FontSize', 16)

%画真值
h30 = scatter(ax3, NaN, NaN, ...                % 数据点的x和y坐标
    50, ...                      % 数据点大小
    [0.5, 0.2, 0.6], ...          % 数据点填充颜色（自定义颜色）
    'filled', ...                 % 填充数据点
    'MarkerEdgeColor', [0.5, 0.2, 0.6], ...   % 数据点边缘颜色（绿色）
    'MarkerFaceAlpha', 0.2, ...   % 设置数据点填充颜色透明度
    'MarkerEdgeAlpha', 1);        % 设置数据点边缘颜色透明度

%画滤波后的
h31 = scatter(ax3, NaN, NaN, ...                % 数据点的x和y坐标
    50, ...                      % 数据点大小
     [183/255, 034/255, 048/255], ...          % 数据点填充颜色（自定义颜色）
    'filled', ...                 % 填充数据点
    'MarkerEdgeColor',  [183/255, 034/255, 048/255], ...   % 数据点边缘颜色（绿色）
    'MarkerFaceAlpha', 0.4, ...   % 设置数据点填充颜色透明度
    'MarkerEdgeAlpha', 1);        % 设置数据点边缘颜色透明度

axis(ax3, 'equal');
xlim(ax3, [1 64])
ylim(ax3, [1 64])

% 横纵坐标显示
xlabel(ax3, 'angle', 'FontSize', 14);
ylabel(ax3, 'range', 'FontSize', 14);

% 关闭x轴和y轴的刻度值
xticks(ax3, []);
yticks(ax3, []);

% 关闭x轴和y轴的刻度线
ax3.XAxis.TickLength = [0 0];
ax3.YAxis.TickLength = [0 0];

% 保留封闭的box
box(ax3, 'on');
axis(ax3, 'equal');

% 添加图标
% legend(ax4, [h31,h32, h33], { 'Resampling Point','Predicted Trajectory', 'Uncertain Region'})
hold(ax3, 'off');

%% 绘制笛卡尔坐标系
hold(ax4, 'on') % 将绘图保持在同一个图中
title(ax4, 'Reconstructed Trajectory', 'FontSize', 16)
[X1,Y1] = meshgrid(1:1:64,1:1:64);
mesh_P = [X1(:),Y1(:)];
Z = zeros(size(X1));
h40 = plot(ax4, NaN, NaN, 'LineWidth', 10,  'Color', [183/255, 034/255, 048/255, 0.4]);        % 设置数据点边缘颜色透明度

% 横纵坐标显示
xlabel(ax4, 'x', 'FontSize', 14);
ylabel(ax4, 'y', 'FontSize', 14);

% 关闭x轴和y轴的刻度值
xticks(ax4, []);
yticks(ax4, []);


% 关闭x轴和y轴的刻度线
ax4.XAxis.TickLength = [0 0];
ax4.YAxis.TickLength = [0 0];

axis(ax4, 'equal');
xlim(ax4, [1 64])
ylim(ax4, [1 64])


% 保留封闭的box
box(ax4, 'on');
hold(ax4, 'off');

%% 绘制AR特征图
hold(ax5, 'on') % 将绘图保持在同一个图中
title(ax5, 'Range-Angle Feature', 'FontSize', 16)

h_img = image(ax5,NaN);
h_scatter = scatter(ax5,NaN, NaN,'MarkerEdgeColor', [0.5 0 0.5],'LineWidth',0.1,'MarkerEdgeAlpha',0.3);
h_center = scatter(ax5,NaN, NaN, 'filled', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'w');
axis(ax5, 'equal');
xlim(ax5, [1 64])
ylim(ax5, [1 64])

% 横纵坐标显示
xlabel(ax5, 'angle', 'FontSize', 14);
ylabel(ax5, 'range', 'FontSize', 14);

% 关闭x轴和y轴的刻度值
xticks(ax5, []);
yticks(ax5, []);

% 关闭x轴和y轴的刻度线
ax5.XAxis.TickLength = [0 0];
ax5.YAxis.TickLength = [0 0];
% 保留封闭的box
box(ax5, 'on');
hold off

btn = uicontrol('Style', 'pushbutton', 'String', 'start', ...
                'Position', [80, 10, 60, 20], 'Callback', @break_blocking_callback);            
 
% 阻塞状态
while true
    % 当用户按下按钮时，回调函数会改变按钮的UserData属性
    if isequal(get(btn, 'UserData'), true)
        break;
    end
    pause(0.5); % 避免过度占用计算资源
end
h = [];    
for i=1:frame_num
    rai_abs = squeeze(range_azimuth(i,:,:));

    set(h_img, 'CData', rai_abs/max(max(rai_abs))*255);
    
    % # 聚类轨迹点
    if scaling_factor(i)>0
        
        [range, angle] = find(rai_abs>=(0.55*max(max(rai_abs))));%找到大于0.55*Max的点
        data = [range, angle];
        [cluster_max_center, cluster_max] = grow_clusters(data);

        %将第一个点当作输入
        tracklen = size(queue_track);
        tracklen = tracklen(1);
        if tracklen==1
            queue_track = cluster_max_center;
        end

        filter_point = kalmanandsmooth_filter(cluster_max_center, queue_track(end,:));
        queue_track = [queue_track;filter_point'];

        %更新图1
        set(h10, 'XData', [get(h10,'XData') i+1], 'YData', [get(h10,'YData') cluster_max_center(2)]);
        set(h12, 'XData', [get(h12,'XData') i+1], 'YData', [get(h12,'YData') queue_track(end,2)]);

        %更新图2
        set(h20, 'XData', [get(h20,'XData') i+1], 'YData', [get(h20,'YData') cluster_max_center(1)]);
        set(h22, 'XData', [get(h22,'XData') i+1], 'YData', [get(h22,'YData') queue_track(end,1)]);


        %更新图5
        set(h_scatter, 'XData',cluster_max(:,2),'YData',cluster_max(:,1));
        set(h_center,  'XData',cluster_max_center(:,2),'YData',cluster_max_center(:,1));
        %更新图3
        set(h30, 'XData', [get(h30,'XData') cluster_max_center(2)], 'YData', [get(h30,'YData') cluster_max_center(1)]);
        set(h31, 'XData', queue_track(:,2), 'YData', queue_track(:,1));
        %更新图4
        queue_track1 = [(queue_track(:, 1)+10).*sin(asin(2*queue_track(:, 2)/64-1)),(queue_track(:, 1)+10).*cos(asin(2*queue_track(:, 2)/64-1))];
        set(h40, 'XData', queue_track1(:,1)+32, 'YData', queue_track1(:,2));

    else
        set(h_scatter, 'XData',NaN,'YData',NaN);
        set(h_center, 'XData',NaN,'YData',NaN);

    end
    
        % 创建一个热图颜色映射
    cmap = colormap(parula);

    % 将图像保存为文件，并使用热图颜色映射
    filename = sprintf('./fig/M_p/frameM_%d.png', i); % 拼接文件名
    filename1 = sprintf('./fig/M_p/frameM1_%d.png', i); % 拼接文件名
%     fig = gcf;  % 获取当前图形句柄
%     axall = findobj(fig, 'Type', 'axes', 'Tag', '');  % 找到所有未标记的轴
%     % 创建一个新的空白图形
%     new_fig = figure('visible', 'off'); % 可以通过将 'off' 改为 'on' 显示新图形
%     % 复制所需子图（在这里是 axall(2)）到新图形
%     new_ax = copyobj(axall(2), new_fig);
%     % 调整新图形的大小和位置
%     set(new_ax, 'Position', [0.1, 0.1, 0.8, 0.8]);
%     saveas(new_ax, filename, 'svg');  % 将图形保存为 SVG 文件
%     close(new_fig);
    
    % 将每帧子图的内容存储到数组中
    image = frame2im(getframe(ax2));
    image = image(2:end,:,:);
    imwrite(image, cmap, filename);
    
    % 将每帧子图的内容存储到数组中
    image = frame2im(getframe(ax1));
    image = image(2:end,:,:);
    imwrite(image, cmap, filename1);
    
%     pause(0.05)
end

%%
% ellipse输出椭圆
% z 是分布
function [ellipse, z, mux, muy, cov, r] = GPR_xy_trajectory(gpr, P, track_xy, section_len)

    len = length(track_xy(:,1));
    
    if len<section_len
        train_X  = [-section_len+len+1:len]'; %#ok<NBRAK>
        
        train_y  = [ones(1,section_len-len)*track_xy(1,1) track_xy(1:len,1)']';%用第一个数填充
        train_y2 = [ones(1,section_len-len)*track_xy(1,2) track_xy(1:len,2)']';
    else
        train_X  = [-section_len+len+1:len]'; %#ok<NBRAK>
        
        train_y  =  track_xy(-section_len+len+1:len,1);
        train_y2 = track_xy(-section_len+len+1:len,2);
    end
    
    train_xy = zeros(2*section_len,1);
    % 交叉重组
    train_xy(1:2:end,1) = train_y;
    train_xy(2:2:end,1) = train_y2;
    
    r = zeros(0, 0);
    % 进行周围采样
    n = 10;
    for j = 1 : length(train_y)
        r((j-1)*n+1:j*n,:) = mvnrnd([train_y(j), train_y2(j)],[1 0;0 1],n);% 设置噪声
    end
    s = invcov_matrix(r(:,1)',r(:,2)');
    
%     s = [1 0;0 1];
%     for j = 1 : length(train_y)
%         s = s + [train_y(j) - mean(train_y),train_y2(j) - mean(train_y2)]'*[train_y(j) - mean(train_y),train_y2(j) - mean(train_y2)];
% %         s = s + u(j);
%     end
   
    
    
    
    test_X = -section_len+len-0.1:0.1:len+1;
    test_X = test_X';
    
    gpr = gpr.fit(train_X, train_xy, s);
    [mu, cov] = gpr.predict(test_X);
    
    mux = mu(1:2:end,:); 
    muy = mu(2:2:end,:); 

%     [ellipse, z] = get_normal_route(s, P, [mux(end) muy(end)]);
    [ellipse, z] = get_normal_route(cov(end-1:end,:), P, [mux(end) muy(end)]);
   
end

function y = kalmanandsmooth_filter(meansure_point, predict_point)
    global P
    global Q
    global R

    P =  P  + Q;         % 卡尔曼公式2
    K = P * inv(P + R); % 卡尔曼公式3
    kf_point = predict_point' + K * (meansure_point'-predict_point');    % 卡尔曼公式4
    P = (eye(2)-K) * P;                       % 卡尔曼公式5
    y = kf_point;
end


function [ellipse, z] = get_normal_route(s, P, next_coordinate)

    [ellipse, angle] = confidence_interval([next_coordinate(1), next_coordinate(2)], s, 0.05);
    L = [P(:,1)-next_coordinate(1),P(:,2)-next_coordinate(2)]*s.*[P(:,1)-next_coordinate(1),P(:,2)-next_coordinate(2)];
    L = sum(L,2);
    len = sqrt(length(L));
    L = reshape(L,[len,len]);
    z = exp(-2.400*L)';%前面的系数是缩放窗口的大小
%     xm = 0:1:63;
%     ym = 0:1:63;
%     ym =ym';
%     xm1 = (xm-next_coordinate(2)).*cos(angle)+(ym-next_coordinate(1)).*sin(angle);
%     ym1 = (ym-next_coordinate(1)).*cos(angle)-(xm-next_coordinate(2)).*sin(angle);
% 
%     x_3db = 0.10*64;
%     y_3db = 0.06*64;
%     n = 3;
%     z = sqrt(1./((xm1/x_3db).^(2*n)+(ym1/y_3db).^(2*n)+1));
end

function s = invcov_matrix(x, y)
    s = pinv(cov([x;y]'));
    if all(s(:)==0)%什么时候出现全零的情况
        s = eye(2);
    end
%     s = eye(2);
end

function [ellipse,angle] = confidence_interval(mu, Sigma, alpha)

    % 计算95%置信区间对应的卡方分布阈值
    chi2_val = chi2inv(1 - alpha/2, 2);
    % 计算置信椭圆的参数
    [V, D] = eig(Sigma);
    [lamda, idx] = sort(diag(D),'descend');
    V = V(:,idx); 
%     V = abs(V);
    
    a = sqrt(chi2_val * abs(lamda(1))); % 长轴
    b = sqrt(chi2_val * abs(lamda(2))); % 短轴
    
    angle = atan2(V(2,1), V(1,1)); % 旋转角度

    % 生成角度向量
    theta = linspace(0, 2 * pi, 100);

    % 计算椭圆上的点
    ellipse_x = 1/a * cos(theta);
    ellipse_y = 1/b * sin(theta);

    % 旋转并平移椭圆
    R = [cos(angle), -sin(angle); sin(angle), cos(angle)];
    ellipse = R * [ellipse_x; ellipse_y];
    ellipse(1, :) = ellipse(1, :) + mu(1);
    ellipse(2, :) = ellipse(2, :) + mu(2);
end


% 回调函数
function break_blocking_callback(src, ~)
    % 将按钮的UserData属性设置为true，以便在主循环中检测到按下按钮的事件
    set(src, 'UserData', true);
end
