%%% Plot firing pattern and spike train
%%% Calculate ISI and CV under different phase
%%% Author: Ye Wang

close all
clear
clc


Position = [1 31.4000 1536 765.6000];

today_date = datestr(now,'yymmdd');

main_dir = 'D:\Memory_integration\';

if ~exist([main_dir,'data\'],'dir')==1
   mkdir([main_dir,'data\']);
end

Save_data_dir = [main_dir,'data\'];

coding_thre = 10;
exc_exc_conn_rate = 0.2;
single_pattern_num = 6;
single_pattern_show_time_new = 10; 
single_pattern_show_time_old = 0;
asso_pattern_show_time_new = 10;
asso_pattern_show_time_old = 0;
validation_show_time = asso_pattern_show_time_new;
pattern_show_time = single_pattern_show_time_new; 
inp_neuron_num = 400;
exc_neuron_num = 800; 

inp_exc_ini = 30;
exc_exc_ini = 5; 
exc_inh_ini = 5;
inh_exc_ini = 25; 
inh_inh_ini = 5;
exc_ref_mu = 20;

pattern_time = 1000;
interval_time = 1000;

A2_plus_exc_novel = 0.4; 
A2_minus_exc_novel = 0.02; 
A3_plus_exc_novel = 0.4;
A3_minus_exc_novel = 0.02; 

trail = 7; % Select a representative trail

stim_file = ['stim\AI_shuffle_show_inp=',num2str(inp_neuron_num),'_pattern_num=',num2str(single_pattern_num),'_single_reps=',num2str(single_pattern_show_time_new),'_asso_reps=',num2str(asso_pattern_show_time_new)...
    '_test_reps=',num2str(validation_show_time),'_pattern_time=',num2str(pattern_time),'_interval_time=',num2str(interval_time),'_trail=',num2str(trail)']; 

stim_path = [main_dir,stim_file,'.mat'];
load(stim_path); 

exc_fr_path = [Save_data_dir,'AI_inp=',num2str(inp_neuron_num),'_',num2str(single_pattern_num),'P_',num2str(single_pattern_show_time_old+single_pattern_show_time_new),'S_',...
    num2str(asso_pattern_show_time_old + asso_pattern_show_time_new),'A_exc_fr_inp-E=',num2str(inp_exc_ini),...
    '_E-E=',num2str(exc_exc_ini),'_E-I=',num2str(exc_inh_ini),'_I-E=',num2str(inh_exc_ini),'_I-I=',num2str(inh_inh_ini),'_pattern_time=',num2str(pattern_time),'_interval_time=',num2str(interval_time),...
    '_LTP_lr=',num2str(A2_plus_exc_novel),'_LTD_lr=',num2str(A2_minus_exc_novel),'_exc_ref=',num2str(exc_ref_mu),'_trail=',num2str(trail),'.mat'];
exc_spike_time_path = [Save_data_dir,'AI_inp=',num2str(inp_neuron_num),'_',num2str(single_pattern_num),'P_',num2str(single_pattern_show_time_old+single_pattern_show_time_new),'S_',...
    num2str(asso_pattern_show_time_old + asso_pattern_show_time_new),'A_exc_spike_time_inp-E=',num2str(inp_exc_ini),...
    '_E-E=',num2str(exc_exc_ini),'_E-I=',num2str(exc_inh_ini),'_I-E=',num2str(inh_exc_ini),'_I-I=',num2str(inh_inh_ini),'_pattern_time=',num2str(pattern_time),'_interval_time=',num2str(interval_time),...
    '_LTP_lr=',num2str(A2_plus_exc_novel),'_LTD_lr=',num2str(A2_minus_exc_novel),'_exc_ref=',num2str(exc_ref_mu),'_trail=',num2str(trail),'.mat'];


load(exc_fr_path);
load(exc_spike_time_path);

pos_S1 = find(strcmp(stim_list, 'S1'));
pos_S2 = find(strcmp(stim_list, 'S2'));
pos_S3 = find(strcmp(stim_list, 'S3'));
pos_S4 = find(strcmp(stim_list, 'S4'));
pos_S5 = find(strcmp(stim_list, 'S5'));
pos_S6 = find(strcmp(stim_list, 'S6'));

pos_A1 = find(strcmp(stim_list, 'A1'));
pos_A2 = find(strcmp(stim_list, 'A2'));
pos_A3 = find(strcmp(stim_list, 'A3'));
pos_A4 = find(strcmp(stim_list, 'A4'));

pos_encoding_end = pos_S6(5);
pos_encoding_start = pos_encoding_end - 1;

pos_vali_end = pos_S1(14);
pos_vali_start = pos_vali_end - 1;

pos_asso_end = pos_A1(4);
pos_asso_start = pos_asso_end - 1;

color_list = [112,193,71;
    237,125,49;
    91,155,213]/255;
hist_bin_num = 10;

% selected encoding phase

k = 1; 
m = 1;

for i = 1:exc_neuron_num
    spike_index = find(exc_spike_time(i,:) > 1000*pos_encoding_start & exc_spike_time(i,:) <= 1000*pos_encoding_end);
    if length(spike_index) > 2 
        for j = 1:length(spike_index)-1
            ISI_encoding(k) = exc_spike_time(i,spike_index(j)+1) - exc_spike_time(i,spike_index(j));
            k = k+1;
            ISI_temp(j) = exc_spike_time(i,spike_index(j)+1) - exc_spike_time(i,spike_index(j));
        end
        CV_ISI_encoding(m) = std(ISI_temp)/mean(ISI_temp);
        m = m + 1;
        clear ISI_temp
    end
end

figure,
set(gcf, 'Position', [680 558 560 204])
subplot(121)
h1 = histogram(ISI_encoding,hist_bin_num);
h1.FaceColor = color_list(1,:);
set(h1,'facealpha',1) 
set(gca,'yscale','log')
xlabel('ISI (ms)')
box off
subplot(122)
h2 = histogram(CV_ISI_encoding,hist_bin_num);
h2.FaceColor = color_list(1,:);
set(h2,'facealpha',1) 
hold on
plot(mean(CV_ISI_encoding),30,'color','k','Marker','v','Markerfacecolor','k','markeredgecolor','k')
xlabel('CV (ISI)')
box off

% selected association phase

k = 1; 
m = 1;

for i = 1:exc_neuron_num
    spike_index = find(exc_spike_time(i,:) > 1000*pos_asso_start & exc_spike_time(i,:) <= 1000*pos_asso_end);
    if length(spike_index) > 2
        for j = 1:length(spike_index)-1
            ISI_asso(k) = exc_spike_time(i,spike_index(j)+1) - exc_spike_time(i,spike_index(j));
            k = k+1;
            ISI_temp(j) = exc_spike_time(i,spike_index(j)+1) - exc_spike_time(i,spike_index(j));
        end
        CV_ISI_asso(m) = std(ISI_temp)/mean(ISI_temp);
        m = m + 1;
        clear ISI_temp
    end
end

figure,
set(gcf, 'Position', [680 558 560 204])
subplot(121)
h1 = histogram(ISI_asso,hist_bin_num);
h1.FaceColor = color_list(2,:);
set(gca,'yscale','log')
set(h1,'facealpha',1) 
xlabel('ISI (ms)')
box off
subplot(122)
h2 = histogram(CV_ISI_asso,hist_bin_num);
h2.FaceColor = color_list(2,:);
set(h2,'facealpha',1) 
hold on
plot(mean(CV_ISI_asso),60,'color','k','Marker','v','Markerfacecolor','k','markeredgecolor','k')
xlabel('CV (ISI)')
box off

% selected validation phase

k = 1; 
m = 1;

for i = 1:exc_neuron_num
    spike_index = find(exc_spike_time(i,:) > 1000*pos_vali_start & exc_spike_time(i,:) <= 1000*pos_vali_end);
    if length(spike_index) > 2
        for j = 1:length(spike_index)-1
            ISI_vali(k) = exc_spike_time(i,spike_index(j)+1) - exc_spike_time(i,spike_index(j));
            k = k+1;
            ISI_temp(j) = exc_spike_time(i,spike_index(j)+1) - exc_spike_time(i,spike_index(j));
        end
        CV_ISI_vali(m) = std(ISI_temp)/mean(ISI_temp);
        m = m + 1;
        clear ISI_temp
    end
end

figure,
set(gcf, 'Position', [680 558 560 204])
subplot(121)
h1 = histogram(ISI_vali,hist_bin_num);
h1.FaceColor = color_list(3,:);
set(h1,'facealpha',1) 
set(gca,'yscale','log')
xlabel('ISI (ms)')
box off
subplot(122)
h2 = histogram(CV_ISI_vali,hist_bin_num);
h2.FaceColor = color_list(3,:);
set(h2,'facealpha',1) 
hold on
plot(mean(CV_ISI_vali),30,'color','k','Marker','v','Markerfacecolor','k','markeredgecolor','k')
xlabel('CV (ISI)')
box off
