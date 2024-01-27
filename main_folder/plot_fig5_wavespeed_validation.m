clear all; close all; clc;
global inputs IC u_data
plot_skip = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Heterogeneous population: D1 = 10D2, A2 = 0.4
% 2. Heterogeneous population: D1 = 10D2, A2 = 0
% 3. Heterogeneous growth population: D1 = D2, A2 = 0
% 4. Homogeneous population
choose_vals = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c_top_sim = [];
c_bot_sim = [];
c_top_data = [];
c_bot_data = [];

save_dir = 'TH_metrics/';
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
save_file = strcat(save_dir,'/wave_speed_sim');
save_file = strcat(save_file,'.mat');
remove_rep = 111:118;
for i=1
    if ~ismember(i,remove_rep)
        data_file = strcat('parameters/full_model/params_cell',num2str(i));
        data_file = strcat(data_file,'.mat');
        % global inputs IC u_data
        load(data_file)

        cell_line_unique = unique(cell_lines,'stable');
        group = find(cell_line_unique == cell_lines(i));
        %%
        delta_t = 0.01;
        t = min(t):delta_t:5;
        x = 0:0.1:5;
        x = x';
        inputs.x = x;
        inputs.t = t;

        minT_idx = floor(0.8*length(t));

        T_interval = t(minT_idx:end);

        u_fmin = u_sum;
        u_data = u_data;
        params_direct = params_direct;
        params_fmin = params_fmin;
        IC = 0.4*exp(-x.^2/(0.1));

        [X,T] = meshgrid(x,t);
        X = X(:,1:plot_skip:end);
        T = T(:,1:plot_skip:end);
        U = u_data(:,1:plot_skip:end);
        X = reshape(X,[],1);
        T = reshape(T,[],1);
        U = reshape(U,[],1);
        
        %%%%%% Test
        %%% Theory fig
        % params_fmin(3) = 5;
        % params_fmin(4) = 2;
        % params_fmin(5) = 0.4;
        % params_fmin(6) = 0.2;
        % params_fmin(7) = 0.6;
        % 
        % params_fmin(1) = 0.1;
        % params_fmin(2) = 2*params_fmin(1);
        % params_fmin(8) = 1.5;
        %%% 1. Heterogeneous population: D1 = 10D2, A2 = 0.4
        if choose_vals == 1
            params_fmin(3) = 2.5;
            params_fmin(4) = 1.5;
            params_fmin(5) = 0.65;
            params_fmin(6) = 0.4;
            params_fmin(7) = 0.5;
            % D1, D2, A2
            params_fmin(1) = 0.007;
            params_fmin(2) = 10*params_fmin(1);
            params_fmin(8) = 0.4;
        %%% 2. Heterogeneous population: D1 = 10D2, A2 = 0
        elseif choose_vals == 2
            params_fmin(3) = 2.5;
            params_fmin(4) = 1.5;
            params_fmin(5) = 0.65;
            params_fmin(6) = 0.4;
            params_fmin(7) = 0.5;
            % D1, D2, A2
            params_fmin(1) = 0.007;
            params_fmin(2) = 10*params_fmin(1);
            params_fmin(8) = 0;
        %%% 3. Heterogeneous growth population: D1 = D2, A2 = 0
        elseif choose_vals == 3
            params_fmin(3) = 2.5;
            params_fmin(4) = 1.5;
            params_fmin(5) = 0.65;
            params_fmin(6) = 0.4;
            params_fmin(7) = 0.5;
            % D1, D2, A2
            params_fmin(1) = 0.007;
            params_fmin(2) = params_fmin(1);
            params_fmin(8) = 0;
        else
        %%% 4. Homogeneous
            params_fmin(3) = 2.5;
            params_fmin(4) = 2.5;
            params_fmin(5) = 0.65;
            params_fmin(6) = 0.65;
            params_fmin(7) = 1;
            
            params_fmin(1) = 0.007;
            params_fmin(2) = params_fmin(1);
            params_fmin(8) = 0;%0.4;
        end
        %%%%%%

        threshold_top = 0.4;%top;
        threshold_bot = 0.1;%bot;
        
        t0_length_str = strcat("t = ",num2str(round(t(floor(0.80*length(t))),2)), " wk");
        t1_length_str = strcat("t = ",num2str(round(t(end),2)), " wk");


        [u_fmin,u1_fmin,u2_fmin] = run_model(params_fmin);
        figure(1); clf;
        hold on
        plot(x,u_fmin(floor(0.80*length(t)),:),'Linewidth',2,'Color',[0 0.4470 0.7410])
        plot(x,u_fmin(end,:),'Linewidth',2,'Color',[0.4660 0.6740 0.1880])
        box on
        yline(threshold_top,'--','Linewidth',1)
        yline(threshold_bot,'--','Linewidth',1)
        set(gca,'FontSize',18)
        xlabel('Radius (mm)')
        ylabel("Density")
        legend(t0_length_str,t1_length_str)
        %%% 1. Heterogeneous population: D1 = 10D2, A2 = 0.4
        if choose_vals == 1
            title('Heterogeneous population: D_1 = 10D_2, A_2 = 0.4')
        %%% 2. Heterogeneous population: D1 = 10D2, A2 = 0
        elseif choose_vals == 2
            title('Heterogeneous population: D_1 = 10D_2, A_2 = 0')
        %%% 3. Heterogeneous growth population: D1 = D2, A2 = 0
        elseif choose_vals == 3
            title('Heterogeneous population: D_1 = D_2, A_2 = 0')
        else
        %%% 4. Homogeneous
            title('Homogeneous population')
        end
        
        xSimTop_tracker = [];
        xSimBot_tracker = [];
        xDataTop_tracker = [];
        xDataBot_tracker = [];
        %%%% Find X
        x_interp = linspace(min(x),max(x),1000);
        for it=minT_idx:length(t)
            %     u_fmin_interp = splineu_fmin(it,:);
            sub_top = abs(u_fmin(it,:) - threshold_top);
            idx_top = find(sub_top == min(sub_top),2);
            xSimTop_tracker = [xSimTop_tracker, x(idx_top)];

            sub_bot = abs(u_fmin(it,:) - threshold_bot);
            idx_bot = find(sub_bot == min(sub_bot),2);
            xSimBot_tracker   = [xSimBot_tracker, x(idx_bot)];

        end

        mSimTop = polyfit(T_interval,xSimTop_tracker,1);
        mSimBot = polyfit(T_interval,xSimBot_tracker,1);

        c_top_sim = [c_top_sim, mSimTop(1)];
        c_bot_sim = [c_bot_sim, mSimBot(1)];
    end

end


function [u_sum,u1,u2] = run_model(params)
global inputs IC model

D_params = [params(1), params(2)];
r_params = [params(3), params(4)];
K_params = [params(5), params(6)];
weights = [params(7), 1-params(7)];
A_params = params(8);
u0 = [weights(1)*IC, weights(2)*IC];


% inputs
x = inputs.x;
t = inputs.t;

try
    [~, uOut] = ode15s(@(t,u)adv_dif_spherical_RHS(x, D_params, r_params, K_params, A_params, u),t,u0);

    numT = length(t);
    numX = length(x);
    numD = length(D_params);
    if numD == 2
        u_temp = reshape(uOut, [numT, numX, numD]);
        u1 = u_temp(:,:,1);
        u2 = u_temp(:,:,2);
        u_sum = u1 + u2;
    else
        u_sum = uOut;
    end
catch
    u_sum = 1e6*ones(length(t),length(x));
end
end
