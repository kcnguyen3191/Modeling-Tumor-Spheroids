clear all; close all; clc;
global inputs IC u_data
plot_skip = 5;
%%%%%%%%%%%% Cell line groups %%%%%%%%%%%%
% Cell line 3013: 1,2,3,4,5,6,7,8 ---- 1:8
% Cell line 3123: 9,10,11,12 ---- 9:12
% Cell line 3118: 13,14,15,16,17,18,19 ---- 13:19
% Cell line 3180: 20,21,22,23,24,25,26,27 ---- 20:27
% Cell line 3289: 28,29,30,31,32,33,34,35 ---- 28:35
% Cell line 3051: 36,37,38,39,40,41,42,43 ---- 36:43
% Cell line 3054: 44,45,46,47,48,49,50,51 ---- 44:51
% Cell line 3279: 52,53,54,55,56,57,58,59 ---- 52:59
% Cell line 3031: 60,61,62,63,64,65,66,67 ---- 60:67
% Cell line 3167: 68,69,70,71,72,73,74,75 ---- 68:75
% Cell line 3291: 76,77,78,79,80,81,82,83 ---- 76:83
% Cell line 3117: 84,85,86,87,88,89       ---- 84:89
% Cell line 3110: 90,91,92,93,94,95,96    ---- 90:96
% Cell line 3021: 97,98,98,99,100,101,102 ---- 97:102
% Cell line 3028: 103,104,105,106,107,108,109,110 ---- 103:110
% Cell line 3179: 111,112,113,114,115,116,117,118 ---- 111:118
% Cell line 3230: 119,120,121,122,123,124,125,126 ---- 119:126
% Cell line 3275: 127,128,129,130,131,132,133,134 ---- 127:134
% Cell line 3086: 135,136,137,138,139,140,141,142 ---- 135:142
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c_top_sim = [];
c_bot_sim = [];
c_top_data = [];
c_bot_data = [];
waveshape = [];
saved_cell_line = [];

% Create save file name
save_dir = 'TH_metrics/';
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
save_file = strcat(save_dir,'/wave_info');
save_file = strcat(save_file,'_sorted.mat');

min_density = 0.02;
max_density = 0.4;
remove_rep = 111:118;
for i=[127:134, 119:126, 135:142, 60:67, 76:83, 68:75, 9:12, 97:102, 90:96, 36:43, 84:89, 44:51, 20:27, 52:59, 103:110, 13:19, 28:35, 1:8] %1:142
    if ~ismember(i,remove_rep)
        %% Load estimated parameter files
        data_file = strcat('parameters/full_model/params_cell',num2str(i));
        data_file = strcat(data_file,'.mat');
        load(data_file)
        
        %% Get unique cell line numbers and find the indidces
        cell_line_unique = unique(cell_lines,'stable');
        group = find(cell_line_unique == cell_lines(i));
        saved_cell_line = [saved_cell_line, cell_lines(i)];
        %%
        delta_t = 0.1;
        t = min(t):delta_t:max(t);
        % x = min(x):0.0001:max(x);
        % x = x';
        inputs.x = x;
        inputs.t = t;
        
        % Choose the time interval to compute the wave speed
        minT_idx = floor(0.75*length(t));
        T_interval = t(minT_idx:end);
        
        % Get fitted results
        u_fmin = u_sum;
        u_data = u_data;
        params_direct = params_direct;
        params_fmin = params_fmin;
        IC = IC;

        [X,T] = meshgrid(x,t);
        X = X(:,1:plot_skip:end);
        T = T(:,1:plot_skip:end);
        U = u_data(:,1:plot_skip:end);
        X = reshape(X,[],1);
        T = reshape(T,[],1);
        U = reshape(U,[],1);

        % [u_fmin,u1_fmin,u2_fmin] = run_model(params_fmin);

        %     threshold_top = top;
        %     threshold_bot = bot;

        figure(1); clf;
        hold on
        plot(x,u_fmin(floor(0.80*length(t)),:),'Linewidth',2,'Color',[0 0.4470 0.7410])
        plot(x,u_fmin(floor(0.90*length(t)),:),'Linewidth',2,'Color',[0.9290 0.6940 0.1250])
        plot(x,u_fmin(end,:),'Linewidth',2,'Color',[0.4660 0.6740 0.1880])
        box on
        set(gca,'FontSize',18)
        xlabel('Position (mm)')
        ylabel("Density")

        max_threshold = 0.8*max(u_fmin(minT_idx,:));
        thresholds = linspace(min_density, max_threshold, 100);
        %%%% Find x function
        WS_sim = [];
        for threshold=thresholds
            xSim_tracker = [];
            % Track x over time
            for it=minT_idx:length(t)
                sub = abs(u_fmin(it,:) - threshold);
                mink_values = mink(sub,1);
                idc = find(ismember(sub, mink_values));
                idx = idc(1);
                xSim_tracker = [xSim_tracker, x(idx)];
            end
            mSim = polyfit(T_interval, xSim_tracker,1);
            WS_sim = [WS_sim, mSim(1)];
            
        end
        WS_sim(WS_sim <= 0) = [];
        c_top_sim = [c_top_sim; min(WS_sim)];
        c_bot_sim = [c_bot_sim; max(WS_sim)];
        
        % Find densities at max and min wave speeds
        top_threshold_idx = find(WS_sim == min(WS_sim));
        top_threshold = thresholds(top_threshold_idx(1));
        bot_threshold_idx = find(WS_sim == max(WS_sim));
        bot_threshold = thresholds(bot_threshold_idx(1));

        sub = abs(u_fmin(end,:) - top_threshold);
        mink_values = mink(sub,1);
        idc = find(ismember(sub, mink_values));
        idx = idc(1);
        top_x = x(idx);

        sub = abs(u_fmin(end,:) - bot_threshold);
        mink_values = mink(sub,1);
        idc = find(ismember(sub, mink_values));
        idx = idc(1);
        bot_x = x(idx);
        
        % Compute the wave shape
        if (top_threshold-bot_threshold)/(top_x-bot_x) == Inf
            i
        end
            waveshape = [waveshape; (top_threshold-bot_threshold)/(top_x-bot_x)];
    end

end
% Compute the wave difference
wavespeed_min = c_top_sim;
wavespeed_max = c_bot_sim;
wavespeed_dif = wavespeed_max - wavespeed_min;
cell_lines(remove_rep) = [];
cell_lines = cell_lines';
save(save_file,'wavespeed_max','wavespeed_min','wavespeed_dif',"c_bot_sim","c_top_sim",'cell_lines','waveshape','saved_cell_line')

function [u_sum,u1,u2] = run_model(params)
global inputs IC
% inputs
x = inputs.x;
t = inputs.t;

D_params = [params(1), params(2)];
r_params = [params(3), params(4)];
K_params = [params(5), params(6)];
weights = [params(7), 1-params(7)];
A_params = params(8);

u0 = [weights(1)*IC, weights(2)*IC];
[~, uOut] = ode15s(@(t,u)adv_dif_spherical_RHS(x, D_params, r_params, K_params, A_params, u),t,u0);


numT = length(t);
numX = length(x);
numD = length(D_params);

u_temp = reshape(uOut, [numT, numX, numD]);
u1 = u_temp(:,:,1);
u2 = u_temp(:,:,2);

u_sum = u1 + u2;
end