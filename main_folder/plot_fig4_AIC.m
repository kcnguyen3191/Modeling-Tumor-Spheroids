clear all; close all; clc;
global inputs IC u_data model
plot_skip = 5;
patient_start = 1;
patient_end = 142;
models = {'FKPP_model',...
    'FKPP_model_wA',...
    'full_model_noA2',...
    'full_model'};
full_model_idx = length(models) + 1;
compared_model_idx = 1+2;

SSE_table = [];
AIC_table = [];
params = {};
for i=patient_start:patient_end
    params{i} = zeros(8,6);
end
remove_rep = [111:118];
for im = 1:length(models)
    if strcmp(models{im},'full_model')
        model = 'full_model';
    elseif strcmp(models{im},'full_model_noA2')
        model = 'full_model_noA2';
    elseif strcmp(models{im},'FKPP_model')
        model = 'FKPP_model';
    elseif strcmp(models{im},'FKPP_model_wA')
        model = 'FKPP_model_wA';
    else
        disp("Enter wrong model name")
        break;
    end

    % Make save directory
    load_dir = strcat('parameters/',model);
    SSE_line = [];
    AIC_line = [];
    cell_lines = [];

    for i=patient_start:patient_end
        if ~ismember(i,remove_rep)
            % Save file
            load_file = strcat(load_dir,'/params_cell',num2str(i),'.mat');
            if exist(load_file,'file')
                parameter_result = load(load_file);

                [m,n] = size(parameter_result.u_data);
                numpar = length(parameter_result.params_fmin);

                cell_lines = [cell_lines, parameter_result.cell_line];

                SSE_fmin = parameter_result.fval_fmin;
                SSE_line = [SSE_line, SSE_fmin];

                AIC_score = m*n*log(SSE_fmin/(m*n)) + m*n*(1 + log(2*pi)) + 2*(numpar + 1);
                AIC_line = [AIC_line, AIC_score];

                if strcmp(models{im},'full_model')
                    params{i}(1:8,im) = parameter_result.params_fmin;
                elseif strcmp(models{im},'full_model_noA2')
                    params{i}(1:numpar,im) = parameter_result.params_fmin;
                elseif strcmp(models{im},'FKPP_model')
                    params{i}([1,3,5],im) = parameter_result.params_fmin;
                elseif strcmp(models{im},'FKPP_model_wA')
                    params{i}([2,4,6,8],im) = parameter_result.params_fmin;
                else
                    disp("Enter wrong model name")
                    break;
                end

            end
        end
    end

    if im == 1
        SSE_table = [cell_lines; SSE_line];
        AIC_table = [cell_lines; AIC_line];
    else
        SSE_table = [SSE_table; SSE_line];
        AIC_table = [AIC_table; AIC_line];
    end
end

SSE_table = SSE_table';
AIC_table = AIC_table';
avg_AIC_table = [];
avg_SSE_table = [];
% cell_lines(remove_rep) = [];
for cell_line = unique(cell_lines,'stable')
    avg_AIC_table = [avg_AIC_table; mean(AIC_table(AIC_table(:,1) == cell_line,:),1)];
    avg_SSE_table = [avg_SSE_table; mean(SSE_table(SSE_table(:,1) == cell_line,:),1)];
end

avg_AIC_table = [avg_AIC_table avg_AIC_table(:,2)-avg_AIC_table(:,full_model_idx)];
avg_SSE_table = [avg_SSE_table avg_SSE_table(:,2)-avg_SSE_table(:,full_model_idx)];

% Add for sorting then remove 
avg_AIC_table = [avg_AIC_table avg_SSE_table(:,6)];
avg_AIC_table = sortrows(avg_AIC_table,7);
avg_AIC_table(:,7) = [];
avg_SSE_table = sortrows(avg_SSE_table,6);

models = {'RD',...
    'ARD',...
    'RD-RD',...
    'RD-ARD'};

xtick_labels = {};
for i=1:length(avg_AIC_table(:,1))
    cell_line=avg_AIC_table(i,1);

    xtick_labels{1,i} = strcat('U',num2str(cell_line),'MG');
end



figure(1)
subplot(2,1,2)
idc = 1:length(unique(cell_lines,'stable'));
for i=1:length(models)
    plot(idc,avg_AIC_table(:,i+1),'linewidth',2)
    hold on
end
hold off
% legend(models)
set(gca,'FontSize',14)
xticks(1:length(unique(cell_lines,'stable')))
xticklabels(xtick_labels)
xlabel('Cell line ID')
ylabel('Average AIC score')
xlim([0,19])
box on
grid on
title('Akaike Information Criterion')

subplot(2,1,1)
idc = 1:length(unique(cell_lines,'stable'));
for i=1:length(models)
    plot(idc,avg_SSE_table(:,i+1),'linewidth',2)
    hold on
end
hold off
legend(models)
set(gca,'FontSize',14)
xticks(1:length(unique(cell_lines,'stable')))
% xticklabels(avg_SSE_table(:,1))
% xlabel('Cell line ID')
set(gca,'Xticklabel',[])
ylabel('Average SSE')
xlim([0,19])
box on
grid on
title('Sum of Squared Errors')

%% Diff
idc = 1:length(unique(cell_lines,'stable'));

figure(5)
subplot(2,1,1)
plot(idc,avg_SSE_table(:,end),'linewidth',2)
set(gca,'FontSize',14)
xticks(1:length(unique(cell_lines,'stable')))
xticklabels(avg_SSE_table(:,1))
xlim([0,19])
xlabel('Cell line ID')
ylabel('Average SSE Difference')
title('Difference in SSE between ARD and RD-ARD models')
box on
grid on

subplot(2,1,2)
plot(idc,avg_AIC_table(:,end),'linewidth',2)
set(gca,'FontSize',14)
xticks(1:length(unique(cell_lines,'stable')))
xticklabels(avg_AIC_table(:,1))
xlabel('Cell line ID')
ylabel('Average AIC score Difference')
xlim([0,19])
box on
grid on
title('Difference in AIC between ARD and RD-ARD models')

save('err_table.mat','idc','avg_SSE_table','avg_AIC_table')