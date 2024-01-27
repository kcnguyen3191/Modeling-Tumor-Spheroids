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
cell_lines = [];
%% Add single cell line value
for i=1:142
% i
data_file = strcat('data/densities_cell',num2str(i));
data_file = strcat(data_file,'.mat');
load(data_file)

full_par_file = strcat('parameters/full_model/params_cell',num2str(i));
full_par_file = strcat(full_par_file,'.mat');
save(full_par_file,'cell_line','-append')

cell_lines = [cell_lines , cell_line];
end

%% Add all single cell line values
for i=1:142
% i
data_file = strcat('data/densities_cell',num2str(i));
data_file = strcat(data_file,'.mat');
load(data_file)

full_par_file = strcat('parameters/full_model/params_cell',num2str(i));
full_par_file = strcat(full_par_file,'.mat');

% global inputs IC u_data
save(full_par_file,'cell_lines','-append')
end
