clear all; close all; clc;
global inputs IC u_data model
plot_skip = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% models = {'full_model',...
%     'full_model_noA2',...
%     'FKPP_model',...
%     'FKPP_model_wA'};
model = {'full_model'};
remove_rep = [111:118];
for i=1:142
    if ~ismember(i,remove_rep)
        i
        data_file = strcat('parameters/',model{1},'/params_cell',num2str(i));
        data_file = strcat(data_file,'.mat')
        % global inputs IC u_data
        load(data_file)
        
        % idc = find(cell_lines == cell_line);
        % rep_idx = find(idc == i);

        %%
        inputs.x = x;
        inputs.t = t;
        u_fmin = u_sum;
        u_data = u_data;
        params_direct = params_direct;
        params_fmin = params_fmin;
        IC = IC;
        length_T = length(t);
        
        [X,T] = meshgrid(x,t);
        X = X(:,1:plot_skip:end);
        T = T(:,1:plot_skip:end);
        U = u_data(:,1:plot_skip:end);
        X = reshape(X,[],1);
        T = reshape(T,[],1);
        U = reshape(U,[],1);        
        
        % [~,u1_fmin,u2_fmin] = run_model(params_fmin);

        idc = find(cell_lines == cell_line);
        rep_idx = find(idc == i);
        %% Figure strings
        if strcmp(model{1},'full_model')
            title_str = "RD-ARD";
        elseif strcmp(model{1},'full_model_noA2')
            title_str = "RD-RD";
        elseif strcmp(model{1},'FKPP_model')
            title_str = "RD";
        elseif strcmp(model{1},'FKPP_model_wA')
            title_str = "ARD";
        end
        subtitle_str = strcat("Cell line: U", num2str(cell_line), "MG, Replicate: ", num2str(rep_idx));
        t0_length_str = strcat("t = ",num2str(round(t(1),2)), " wk");
        t1_length_str = strcat("t = ",num2str(round(0.25*t(end),2)), " wk");
        t2_length_str = strcat("t = ",num2str(round(0.5*t(end),2)), " wk");
        t3_length_str = strcat("t = ",num2str(round(0.75*t(end),2)), " wk");
        tf_length_str = strcat("t = ",num2str(round(t(end),2)), " wk");
        fitting_surface_fig_path = strcat('figures/fitting_results/',model{1},'/',model{1},'_surface',num2str(i),...
            '_cell_line',num2str(cell_line),'_rep',num2str(rep_idx));
        u1_surface_fig_path = strcat('figures/fitting_results/',model{1},'/',model{1},'_u1_surface',num2str(i),...
            '_cell_line',num2str(cell_line),'_rep',num2str(rep_idx));
        u2_surface_fig_path = strcat('figures/fitting_results/',model{1},'/',model{1},'_u2_surface',num2str(i),...
            '_cell_line',num2str(cell_line),'_rep',num2str(rep_idx));
        fitting_curve_fig_path = strcat('figures/fitting_results/',model{1},'/',model{1},'_curve',num2str(i),...
            '_cell_line',num2str(cell_line),'_rep',num2str(rep_idx));

        %% Plotting
        figure(1); clf;
        mesh(x,t,u_fmin)
        box on
        hold on
        scatter3(X,T,U,'filled','MarkerFaceColor',[0 0 0],'LineWidth',0.5);
        hold off
        xlabel('Position (mm)')
        ylabel('Time (weeks)')
        zlabel('Densitiy')
        set(gca,'FontSize',15)
        title(title_str,'FontSize',22)
        subtitle(subtitle_str,'FontSize',20)
        xlim([0,2])
        zlim([0,0.5])
        view(-3,30)
        saveas(gcf,fitting_surface_fig_path,'epsc')
        
        figure(4); clf;
        hold on
        plot(x,u_fmin(1,:),'Linewidth',2,'Color',[0 0.4470 0.7410])
        plot(x,u_fmin(floor(length_T/4),:),'Linewidth',2,'Color',[0.8500 0.3250 0.0980])
        plot(x,u_fmin(floor(length_T/2),:),'Linewidth',2,'Color',[0.9290 0.6940 0.1250])
        plot(x,u_fmin(floor(3*length_T/4),:),'Linewidth',2,'Color',[0.4940 0.1840 0.5560])
        plot(x,u_fmin(end,:),'Linewidth',2,'Color',[0.4660 0.6740 0.1880])
        box on
        set(gca,'FontSize',18)
        title(title_str,'FontSize',22)
        subtitle(subtitle_str,'FontSize',20)
        xlabel('Position (mm)')
        ylabel("Density")
        plot(x,u_data(1,:),'.-','Color',[0 0.4470 0.7410])
        plot(x,u_data(floor(length_T/4),:),'.-','Color',[0.8500 0.3250 0.0980])
        plot(x,u_data(floor(length_T/2),:),'.-','Color',[0.9290 0.6940 0.1250])
        plot(x,u_data(floor(3*length_T/4),:),'.-','Color',[0.4940 0.1840 0.5560])
        plot(x,u_data(length_T,:),'.-','Color',[0.4660 0.6740 0.1880])
        hold off
        legend(t0_length_str,t1_length_str,t2_length_str,t3_length_str,tf_length_str)
        xlim([0,2])
        ylim([0,0.5])
        saveas(gcf,fitting_curve_fig_path,'epsc')
    end
end

function [u_sum,u1,u2] = run_model(params)
global inputs IC model


if strcmp(model,'full_model')
    D_params = [params(1), params(2)];
    r_params = [params(3), params(4)];
    K_params = [params(5), params(6)];
    weights = [params(7), 1-params(7)];
    A_params = params(8);
    u0 = [weights(1)*IC, weights(2)*IC];
elseif strcmp(model,'full_model_noA2')
    D_params = [params(1), params(2)];
    r_params = [params(3), params(4)];
    K_params = [params(5), params(6)];
    weights = [params(7), 1-params(7)];
    A_params = 0;
    u0 = [weights(1)*IC, weights(2)*IC];
elseif strcmp(model,'FKPP_model')
    D_params = params(1);
    r_params = params(2);
    K_params = params(3);
    A_params = 0;
    u0 = IC;
elseif strcmp(model,'FKPP_model_wA')
    D_params = params(1);
    r_params = params(2);
    K_params = params(3);
    A_params = params(4);
    u0 = IC;
else
    disp("Enter wrong model name")
end


% inputs
x = inputs.x;
t = inputs.t;

try
    [~, uOut] = ode45(@(t,u)adv_dif_spherical_RHS(x, D_params, r_params, K_params, A_params, u),t,u0);

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