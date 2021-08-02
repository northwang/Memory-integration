%%% Plot weight change of selected coding neurons
%%% The selected coding neurons are always coding each pattern across time
%%% Author: Ye Wang

close all
clear
clc

today_date = datestr(now,'yymmdd');

main_dir = 'D:\NIMI\Programs\Cell_assembly\';

if ~exist([main_dir,'data_',num2str(today_date)],'dir')==1
   mkdir([main_dir,'data_',num2str(today_date)]);
end

Save_data_dir = [main_dir,'data_',num2str(today_date),'\'];

coding_thre = 40;
exc_exc_conn_rate = 0.2;
single_pattern_num = 6;
single_pattern_show_time_new = 10;
single_pattern_show_time_old = 0;
asso_pattern_show_time_new = 10;
asso_pattern_show_time_old = 0;
validation_show_time = asso_pattern_show_time_new;
total_show_time = single_pattern_show_time_new + validation_show_time;
inp_neuron_num = 400;
exc_neuron_num = 800; 

inp_exc_ini = 30;
exc_exc_ini = 5; 
exc_inh_ini = 5;
inh_exc_ini = 25; 
inh_inh_ini = 5;
exc_ref_mu = 20;

max_exc_exc_weight = 10*exc_exc_ini;
min_exc_exc_weight = 0.1*exc_exc_ini;

pattern_time = 1000;
interval_time = 1000;

A2_plus_exc_novel = 0.4; 
A2_minus_exc_novel = 0.02; 
A3_plus_exc_novel = 0.4;
A3_minus_exc_novel = 0.02; 

trail = 7;

stim_file = ['stim_210618\AI_shuffle_show_inp=',num2str(inp_neuron_num),'_pattern_num=',num2str(single_pattern_num),'_single_reps=',num2str(single_pattern_show_time_new),'_asso_reps=',num2str(asso_pattern_show_time_new)...
    '_test_reps=',num2str(validation_show_time),'_pattern_time=',num2str(pattern_time),'_interval_time=',num2str(interval_time),'_trail=',num2str(trail)]; 

stim_path = [main_dir,stim_file,'.mat'];
load(stim_path); % ¶ÁÈ¡´Ì¼¤ÎÄ¼þ

exc_fr_path = [Save_data_dir,'AI_inp=',num2str(inp_neuron_num),'_',num2str(single_pattern_num),'P_',num2str(single_pattern_show_time_old+single_pattern_show_time_new),'S_',...
    num2str(asso_pattern_show_time_old + asso_pattern_show_time_new),'A_exc_fr_inp-E=',num2str(inp_exc_ini),...
    '_E-E=',num2str(exc_exc_ini),'_E-I=',num2str(exc_inh_ini),'_I-E=',num2str(inh_exc_ini),'_I-I=',num2str(inh_inh_ini),'_pattern_time=',num2str(pattern_time),'_interval_time=',num2str(interval_time),...
    '_LTP_lr=',num2str(A2_plus_exc_novel),'_LTD_lr=',num2str(A2_minus_exc_novel),'_exc_ref=',num2str(exc_ref_mu),'_trail=',num2str(trail),'.mat'];
EE_weight_path = [Save_data_dir,'AI_inp=',num2str(inp_neuron_num),'_',num2str(single_pattern_num),'P_',num2str(single_pattern_show_time_old+single_pattern_show_time_new),'S_',...
    num2str(asso_pattern_show_time_old + asso_pattern_show_time_new),'A_EE_weight_inp-E=',num2str(inp_exc_ini),...
    '_E-E=',num2str(exc_exc_ini),'_E-I=',num2str(exc_inh_ini),'_I-E=',num2str(inh_exc_ini),'_I-I=',num2str(inh_inh_ini),'_pattern_time=',num2str(pattern_time),'_interval_time=',num2str(interval_time),...
    '_LTP_lr=',num2str(A2_plus_exc_novel),'_LTD_lr=',num2str(A2_minus_exc_novel),'_exc_ref=',num2str(exc_ref_mu),'_trail=',num2str(trail),'.mat'];

load(exc_fr_path);
load(EE_weight_path); 

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

A_sel = A_coding_neuron{1,1};
B_sel = B_coding_neuron{1,1};
C_sel = C_coding_neuron{1,1};
D_sel = D_coding_neuron{1,1};
E_sel = E_coding_neuron{1,1};
F_sel = F_coding_neuron{1,1};
O_sel = O_coding_neuron{1,1};

for i = 1:total_show_time
    A_sel = intersect(A_sel,A_coding_neuron{i,1});
    B_sel = intersect(B_sel,B_coding_neuron{i,1});
    C_sel = intersect(C_sel,C_coding_neuron{i,1});
    D_sel = intersect(D_sel,D_coding_neuron{i,1});
    E_sel = intersect(E_sel,E_coding_neuron{i,1});
    F_sel = intersect(F_sel,F_coding_neuron{i,1});
    O_sel = intersect(O_sel,O_coding_neuron{i,1});
end

for i = 1:size(exc_exc_weight_save,3)
    weight_A2A(i) = mean_weight(exc_exc_weight_save,A_sel,A_sel,i);
    weight_A2B(i) = mean_weight(exc_exc_weight_save,A_sel,B_sel,i);
    weight_A2C(i) = mean_weight(exc_exc_weight_save,A_sel,C_sel,i);
    weight_A2D(i) = mean_weight(exc_exc_weight_save,A_sel,D_sel,i);
    weight_A2O(i) = mean_weight(exc_exc_weight_save,A_sel,O_sel,i);
end

figure,
set(gcf,'Position',[11.4000 487.4000 1528 241.2000])
plot(weight_A2A)
hold on
plot(weight_A2B)
hold on
plot(weight_A2C)
hold on
plot(weight_A2D)
hold on
plot(weight_A2O)
legend('A to A (inner)','A to B (0¡ã)','A to C (1¡ã)','A to D (inf)','A to background (others)')
xlim([0,size(exc_exc_weight_save,3)])
xlabel('Time (s)')
ylabel('Weight')
box off


