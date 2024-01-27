clear all; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TH_folder = 'TH_metrics/';
WS_file = strcat(TH_folder,'wave_info_sorted.mat');
load(WS_file);

load('err_table.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove_rep = 111:118;
% c_bot_sim(remove_rep) = [];
% c_top_sim(remove_rep) = [];
% cell_lines(remove_rep)= [];
cell_line_unique = avg_SSE_table(:,1)';%[3275 3123 3291 3086 3230 3117 3031 3051 3167 3180 3110 3021 3054 3118 3279 3028 3289 3013];%unique(cell_lines,'stable');
group_count = 1;

g = [];
xtick_labels = {};
for i=1:length(cell_line_unique)
    cell_line=cell_line_unique(i);
    rep_idc = find(cell_lines == cell_line);

    cell_idx = find(cell_line_unique == cell_line);
    g = [g; cell_idx*ones(length(rep_idc),1)];

    xtick_labels{1,i} = strcat('U',num2str(cell_line),'MG');
end

scale_x = 1;
figure(1)
subplot(4,1,3)
boxchart(g*scale_x,c_bot_sim-c_top_sim)
hold on
box on
grid on
xticks(scale_x*(1:length(cell_line_unique)))
% xticklabels(xtick_labels)
set(gca,'Xticklabel',[])
ylim([0, 0.5])
xlim([0, 19])
set(gca,'FontSize',14)
title('c_{diff}')
% xlabel('Cell line ID')
ylabel("Wave speed")

subplot(4,1,4)
boxchart(g*scale_x,waveshape)
hold on
box on
grid on
xticks(scale_x*(1:length(cell_line_unique)))
xticklabels(xtick_labels)
% ylim([0, 0.5])
xlim([0, 19])
set(gca,'FontSize',14)
title('c_{shape}')
xlabel('Cell line ID')
ylabel("Wave shape")

idc = 1:length(unique(cell_lines,'stable'));

subplot(4,1,1)
plot(idc,avg_SSE_table(:,end),'linewidth',2)
set(gca,'FontSize',14)
xticks(1:length(unique(cell_lines,'stable')))
% xticklabels(xtick_labels)
set(gca,'Xticklabel',[])
xlim([0,19])
% xlabel('Cell line ID')
ylabel('\Delta SSE')
title('Difference in SSE between ARD and RD-ARD models')
box on
grid on

subplot(4,1,2)
plot(idc,avg_AIC_table(:,end),'linewidth',2)
set(gca,'FontSize',14)
xticks(1:length(unique(cell_lines,'stable')))
% xticklabels(xtick_labels)
set(gca,'Xticklabel',[])
xlim([0,19])
% xlabel('Cell line ID')
ylabel('\Delta AIC')
title('Difference in AIC between ARD and RD-ARD models')
box on
grid on


figure(2)
boxchart(g*scale_x,c_bot_sim)
hold on
box on
grid on
boxchart(scale_x*g,c_top_sim)
xticks(scale_x*(1:length(cell_line_unique)))
xticklabels(cell_line_unique)
set(gca,'FontSize',14)
title('Boxplot: c_{max} and c_{min}')
xlabel('Cell line ID')
ylabel("Wavespeed")
xlim([0, 19])
legend('c_{max}','c_{min}')