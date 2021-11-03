%%% Plot weight change of selected coding neurons
%%% The selected coding neurons are always coding each pattern across time
%%% Author: Ye Wang

close all
clear
clc

main_dir = 'D:\NIMI\Programs\Cell_assembly\';

if ~exist([main_dir,'data\'],'dir')==1
   mkdir([main_dir,'data\']);
end

Save_data_dir = [main_dir,'data\'];

coding_thre = 40;
exc_exc_conn_rate = 0.2;
single_pattern_num = 6;
single_pattern_show_time_new = 10; % 这次学习呈现的次数
single_pattern_show_time_old = 0;
asso_pattern_show_time_new = 10;
asso_pattern_show_time_old = 0;
validation_show_time = asso_pattern_show_time_new;
total_show_time = single_pattern_show_time_new + validation_show_time;
inp_neuron_num = 400;
exc_neuron_num = 800; 

inp_exc_ini = 30;
exc_exc_ini = 10; % 10
exc_inh_ini = 5;
inh_exc_ini = 30; % 30
inh_inh_ini = 5;
exc_ref_mu = 20;

max_exc_exc_weight = 10*exc_exc_ini;
min_exc_exc_weight = 0.1*exc_exc_ini;

pattern_time = 1000;
interval_time = 1000;
Nskip = 1000; % how often (in number of timesteps) to save weight

A2_plus_exc_novel = 0.4; 
A2_minus_exc_novel = 0.02; 
A3_plus_exc_novel = 0.4;
A3_minus_exc_novel = 0.02; 

encoding_time_start = 109;
encoding_time_end = 121;

vali_time_start = 301;
vali_time_end = 320;

trial = 6；

stim_file = ['stim_211009\AI_shuffle_show_inp=',num2str(inp_neuron_num),'_pattern_num=',num2str(single_pattern_num),'_single_reps=',num2str(single_pattern_show_time_new),'_asso_reps=',num2str(asso_pattern_show_time_new)...
    '_test_reps=',num2str(validation_show_time),'_pattern_time=',num2str(pattern_time),'_interval_time=',num2str(interval_time),'_trail=',num2str(trial)]; 

stim_path = [main_dir,stim_file,'.mat'];
load(stim_path); % 读取刺激文件

exc_fr_path = [Save_data_dir,'AI_inp=',num2str(inp_neuron_num),'_',num2str(single_pattern_num),'P_',num2str(single_pattern_show_time_old+single_pattern_show_time_new),'S_',...
    num2str(asso_pattern_show_time_old + asso_pattern_show_time_new),'A_exc_fr_inp-E=',num2str(inp_exc_ini),...
    '_E-E=',num2str(exc_exc_ini),'_E-I=',num2str(exc_inh_ini),'_I-E=',num2str(inh_exc_ini),'_I-I=',num2str(inh_inh_ini),'_pattern_time=',num2str(pattern_time),'_interval_time=',num2str(interval_time),...
    '_LTP_lr=',num2str(A2_plus_exc_novel),'_LTD_lr=',num2str(A2_minus_exc_novel),'_exc_ref=',num2str(exc_ref_mu),'_trail=',num2str(trial),'.mat'];
EE_weight_path = [Save_data_dir,'AI_inp=',num2str(inp_neuron_num),'_',num2str(single_pattern_num),'P_',num2str(single_pattern_show_time_old+single_pattern_show_time_new),'S_',...
    num2str(asso_pattern_show_time_old + asso_pattern_show_time_new),'A_EE_weight_inp-E=',num2str(inp_exc_ini),...
    '_E-E=',num2str(exc_exc_ini),'_E-I=',num2str(exc_inh_ini),'_I-E=',num2str(inh_exc_ini),'_I-I=',num2str(inh_inh_ini),'_pattern_time=',num2str(pattern_time),'_interval_time=',num2str(interval_time),...
    '_LTP_lr=',num2str(A2_plus_exc_novel),'_LTD_lr=',num2str(A2_minus_exc_novel),'_exc_ref=',num2str(exc_ref_mu),'_trail=',num2str(trial),'.mat'];

load(exc_fr_path);
load(EE_weight_path); % 读取weight matrix文件
    
% 获得每个刺激出现的时间

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


%% Weight change between patterns    
for i = 1:total_show_time
    A_coding_neuron{i,1} = find(exc_fr(:,pos_S1(i)) > coding_thre);
    B_coding_neuron{i,1} = find(exc_fr(:,pos_S2(i)) > coding_thre);
    C_coding_neuron{i,1} = find(exc_fr(:,pos_S3(i)) > coding_thre);
    D_coding_neuron{i,1} = find(exc_fr(:,pos_S4(i)) > coding_thre);
    E_coding_neuron{i,1} = find(exc_fr(:,pos_S5(i)) > coding_thre);
    F_coding_neuron{i,1} = find(exc_fr(:,pos_S6(i)) > coding_thre);
    
    O_coding_neuron{i,1} = setdiff(1:exc_neuron_num,[A_coding_neuron{i,1};B_coding_neuron{i,1};C_coding_neuron{i,1};
            D_coding_neuron{i,1};E_coding_neuron{i,1};F_coding_neuron{i,1}]); % other neurons
end


for i = 1:total_show_time
    if i <= 1
        j_start = 1;
    else
        j_start = max([pos_S1(i-1),pos_S2(i-1),pos_S3(i-1),pos_S4(i-1),pos_S5(i-1),pos_S6(i-1)]) + 1;
    end
    for j = j_start:max([pos_S1(i),pos_S2(i),pos_S3(i),pos_S4(i),pos_S5(i),pos_S6(i)])   
        weight_A2A(j) = mean_weight(exc_exc_weight_save,A_coding_neuron{i,1},A_coding_neuron{i,1},j);
        weight_A2B(j) = mean_weight(exc_exc_weight_save,A_coding_neuron{i,1},B_coding_neuron{i,1},j);
        weight_A2C(j) = mean_weight(exc_exc_weight_save,A_coding_neuron{i,1},C_coding_neuron{i,1},j);
        weight_A2D(j) = mean_weight(exc_exc_weight_save,A_coding_neuron{i,1},D_coding_neuron{i,1},j);
        weight_A2O(j) = mean_weight(exc_exc_weight_save,A_coding_neuron{i,1},O_coding_neuron{i,1},j);
        
        A_fr(j) = mean(exc_fr(A_coding_neuron{i,1},j));
        B_fr(j) = mean(exc_fr(B_coding_neuron{i,1},j));
        C_fr(j) = mean(exc_fr(C_coding_neuron{i,1},j));
        D_fr(j) = mean(exc_fr(D_coding_neuron{i,1},j));
        E_fr(j) = mean(exc_fr(E_coding_neuron{i,1},j));
        F_fr(j) = mean(exc_fr(F_coding_neuron{i,1},j));
    end
end

zero_inf_ratio(trial) = weight_A2B(end)/weight_A2D(end);
one_inf_ratio(trial) = weight_A2C(end)/weight_A2D(end);

color_list = [0.47,0.67,0.19;
    0.00,0.45,0.74;
    0.85,0.33,0.1;
    0.49,0.18,0.56;
    0.93,0.69,0.13;
    1,0.07,0.65];

figure,
set(gcf,'Position',[11.4000 487.4000 1528 241.2000])
plot(weight_A2A,'Color',color_list(1,:),'LineWidth',2);
hold on
plot(weight_A2B,'Color',color_list(2,:),'LineWidth',2);
hold on
plot(weight_A2C,'Color',color_list(3,:),'LineWidth',2);
hold on
plot(weight_A2D,'Color',color_list(4,:),'LineWidth',2);
hold on
plot(weight_A2O,'Color',color_list(5,:),'LineWidth',2);
legend('A to A (inner)','A to B (0°)','A to C (1°)','A to X (inf)','A to background (others)','FontName','Arial','FontSize',12,'NumColumns',5)
legend('boxoff')
set(gca,'FontName','Arial','FontSize',12)
xlim([0,size(exc_exc_weight_save,3)])
xlabel('Time (s)','FontName','Arial','FontSize',12)
ylabel('Weight','FontName','Arial','FontSize',12)
box off
% title(['trial=',int2str(trial)])
saveas(gcf,[Save_figure_dir,'Weight_change_3F_trial=',num2str(trial),'.jpg']);  %保存当前窗口的图像

figure,
set(gcf,'Position',  [366 321 830 160])
plot(encoding_time_start:encoding_time_end,A_fr(encoding_time_start:encoding_time_end),'Color',color_list(1,:),'LineWidth',2);
hold on
plot(encoding_time_start:encoding_time_end,B_fr(encoding_time_start:encoding_time_end),'Color',color_list(2,:),'LineWidth',2);
hold on
plot(encoding_time_start:encoding_time_end,C_fr(encoding_time_start:encoding_time_end),'Color',color_list(3,:),'LineWidth',2);
hold on
plot(encoding_time_start:encoding_time_end,D_fr(encoding_time_start:encoding_time_end),'Color',color_list(4,:),'LineWidth',2);
hold on
plot(encoding_time_start:encoding_time_end,E_fr(encoding_time_start:encoding_time_end),'Color',color_list(5,:),'LineWidth',2);
hold on
plot(encoding_time_start:encoding_time_end,F_fr(encoding_time_start:encoding_time_end),'Color',color_list(6,:),'LineWidth',2);
legend('A','B','C','X','Y','Z','FontName','Arial','FontSize',12,'NumColumns',6)
legend('boxoff')
set(gca,'FontName','Arial','FontSize',12)
xlim([encoding_time_start,encoding_time_end])
set(gca,'xtick',[],'xticklabel',[],'xcolor','w')
ylim([0,100])
% xlabel('Time (s)','FontName','Arial','FontSize',12)
ylabel('Activity (Hz)','FontName','Arial','FontSize',12)
box off
% title(['trial=',int2str(trial)])
% saveas(gcf,[Save_figure_dir,'Mean_fr_3F_trial=',num2str(trial),'.jpg']);  %保存当前窗口的图像

figure,
set(gcf,'Position',  [306 224 1382 160])
plot(vali_time_start:vali_time_end+1,[A_fr(vali_time_start:vali_time_end),0],'Color',color_list(1,:),'LineWidth',2);
hold on
plot(vali_time_start:vali_time_end+1,[B_fr(vali_time_start:vali_time_end),0],'Color',color_list(2,:),'LineWidth',2);
hold on
plot(vali_time_start:vali_time_end+1,[C_fr(vali_time_start:vali_time_end),0],'Color',color_list(3,:),'LineWidth',2);
hold on
plot(vali_time_start:vali_time_end+1,[D_fr(vali_time_start:vali_time_end),0],'Color',color_list(4,:),'LineWidth',2);
hold on
plot(vali_time_start:vali_time_end+1,[E_fr(vali_time_start:vali_time_end),0],'Color',color_list(5,:),'LineWidth',2);
hold on
plot(vali_time_start:vali_time_end+1,[F_fr(vali_time_start:vali_time_end),0],'Color',color_list(6,:),'LineWidth',2);
legend('A','B','C','X','Y','Z','FontName','Arial','FontSize',12,'NumColumns',6)
legend('boxoff')
set(gca,'FontName','Arial','FontSize',12)
xlim([vali_time_start,vali_time_end+1])
set(gca,'xtick',[],'xticklabel',[],'xcolor','w')
ylim([0,100])
% xlabel('Time (s)','FontName','Arial','FontSize',12)
ylabel('Activity (Hz)','FontName','Arial','FontSize',12)
box off

