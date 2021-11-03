%%% Load the simulation results of each trail and analyze weight and overlap 
%%% Author: Ye Wang

close all
clear
clc

main_dir = 'D:\Memory_integration\';

if ~exist([main_dir,'data\'],'dir')==1
   mkdir([main_dir,'data\']);
end

Save_data_dir = [main_dir,'data\'];

coding_thre = 40;
single_pattern_num = 6;
single_pattern_show_time_new = 10; % 这次学习呈现的次数
single_pattern_show_time_old = 0;
asso_pattern_show_time_new = 10;
asso_pattern_show_time_old = 0;
validation_show_time = asso_pattern_show_time_new;
inp_neuron_num = 400;
exc_neuron_num = 800; 

inp_exc_ini = 30;
exc_exc_ini = 10; 
exc_inh_ini = 5;
inh_exc_ini = 30; 
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

trail_num = 30;

overlap_AB = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
overlap_AC = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
overlap_AD = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
overlap_AE = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
overlap_AF = zeros(trail_num,single_pattern_show_time_new+validation_show_time);

overlap_BC = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
overlap_BD = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
overlap_BE = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
overlap_BF = zeros(trail_num,single_pattern_show_time_new+validation_show_time);

overlap_CD = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
overlap_CE = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
overlap_CF = zeros(trail_num,single_pattern_show_time_new+validation_show_time);

overlap_DE = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
overlap_DF = zeros(trail_num,single_pattern_show_time_new+validation_show_time);

overlap_EF = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
mean_overlap = zeros(trail_num,single_pattern_show_time_new+validation_show_time);

weight_A2A = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_A2B = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_A2C = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_A2D = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_A2E = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_A2F = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_A2O = zeros(trail_num,single_pattern_show_time_new+validation_show_time);

weight_B2A = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_B2B = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_B2C = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_B2D = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_B2E = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_B2F = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_B2O = zeros(trail_num,single_pattern_show_time_new+validation_show_time);

weight_C2A = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_C2B = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_C2C = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_C2D = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_C2E = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_C2F = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_C2O = zeros(trail_num,single_pattern_show_time_new+validation_show_time);

weight_D2A = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_D2B = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_D2C = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_D2D = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_D2E = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_D2F = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_D2O = zeros(trail_num,single_pattern_show_time_new+validation_show_time);

weight_E2A = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_E2B = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_E2C = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_E2D = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_E2E = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_E2F = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_E2O = zeros(trail_num,single_pattern_show_time_new+validation_show_time);

weight_F2A = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_F2B = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_F2C = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_F2D = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_F2E = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_F2F = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_F2O = zeros(trail_num,single_pattern_show_time_new+validation_show_time);

weight_O2A = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_O2B = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_O2C = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_O2D = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_O2E = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_O2F = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
weight_O2O = zeros(trail_num,single_pattern_show_time_new+validation_show_time);

A_coding_neuron_num = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
A_coding_neuron_fr = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
B_coding_neuron_num = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
B_coding_neuron_fr = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
C_coding_neuron_num = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
C_coding_neuron_fr = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
D_coding_neuron_num = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
D_coding_neuron_fr = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
E_coding_neuron_num = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
E_coding_neuron_fr = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
F_coding_neuron_num = zeros(trail_num,single_pattern_show_time_new+validation_show_time);
F_coding_neuron_fr = zeros(trail_num,single_pattern_show_time_new+validation_show_time);

for trail = 1:trail_num
    stim_file = ['stim\AI_shuffle_show_inp=',num2str(inp_neuron_num),'_pattern_num=',num2str(single_pattern_num),'_single_reps=',num2str(single_pattern_show_time_new),'_asso_reps=',num2str(asso_pattern_show_time_new)...
        '_test_reps=',num2str(validation_show_time),'_pattern_time=',num2str(pattern_time),'_interval_time=',num2str(interval_time),'_trail=',num2str(trail)]; 

    stim_path = [main_dir,stim_file,'.mat'];
    load(stim_path);

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

    % Get the show time of each stimulus

    pos_S1 = find(strcmp(stim_list, 'S1'));
    pos_S2 = find(strcmp(stim_list, 'S2'));
    pos_S3 = find(strcmp(stim_list, 'S3'));
    pos_S4 = find(strcmp(stim_list, 'S4'));
    pos_S5 = find(strcmp(stim_list, 'S5'));
    pos_S6 = find(strcmp(stim_list, 'S6'));

    pos_S = [pos_S1,pos_S2,pos_S3,pos_S4,pos_S5,pos_S6];


    %% Weight change between patterns    
    for i = 1:length(pos_S1) 

        A_coding_neuron{trail,i} = find(exc_fr(:,pos_S1(i)) > coding_thre);
        B_coding_neuron{trail,i} = find(exc_fr(:,pos_S2(i)) > coding_thre);
        C_coding_neuron{trail,i} = find(exc_fr(:,pos_S3(i)) > coding_thre);
        D_coding_neuron{trail,i} = find(exc_fr(:,pos_S4(i)) > coding_thre);
        E_coding_neuron{trail,i} = find(exc_fr(:,pos_S5(i)) > coding_thre);
        F_coding_neuron{trail,i} = find(exc_fr(:,pos_S6(i)) > coding_thre);
        
        O_coding_neuron{trail,i} = setdiff(1:exc_neuron_num,[A_coding_neuron{trail,i};B_coding_neuron{trail,i};C_coding_neuron{trail,i};
            D_coding_neuron{trail,i};E_coding_neuron{trail,i};F_coding_neuron{trail,i}]); % other neurons

        A_coding_neuron_num(trail,i) = length(A_coding_neuron{trail,i});
        A_coding_neuron_fr(trail,i) = mean(exc_fr(A_coding_neuron{trail,i},pos_S1(i)));
        B_coding_neuron_num(trail,i) = length(B_coding_neuron{trail,i});
        B_coding_neuron_fr(trail,i) = mean(exc_fr(B_coding_neuron{trail,i},pos_S2(i)));
        C_coding_neuron_num(trail,i) = length(C_coding_neuron{trail,i});
        C_coding_neuron_fr(trail,i) = mean(exc_fr(C_coding_neuron{trail,i},pos_S3(i)));
        D_coding_neuron_num(trail,i) = length(D_coding_neuron{trail,i});
        D_coding_neuron_fr(trail,i) = mean(exc_fr(D_coding_neuron{trail,i},pos_S4(i)));
        E_coding_neuron_num(trail,i) = length(E_coding_neuron{trail,i});
        E_coding_neuron_fr(trail,i) = mean(exc_fr(E_coding_neuron{trail,i},pos_S5(i)));
        F_coding_neuron_num(trail,i) = length(F_coding_neuron{trail,i});
        F_coding_neuron_fr(trail,i) = mean(exc_fr(F_coding_neuron{trail,i},pos_S6(i)));

        overlap_AB(trail,i) = coding_overlap(A_coding_neuron{trail,i}, B_coding_neuron{trail,i});
        overlap_AC(trail,i) = coding_overlap(A_coding_neuron{trail,i}, C_coding_neuron{trail,i});
        overlap_AD(trail,i) = coding_overlap(A_coding_neuron{trail,i}, D_coding_neuron{trail,i});
        overlap_AE(trail,i) = coding_overlap(A_coding_neuron{trail,i}, E_coding_neuron{trail,i});
        overlap_AF(trail,i) = coding_overlap(A_coding_neuron{trail,i}, F_coding_neuron{trail,i});

        overlap_BC(trail,i) = coding_overlap(B_coding_neuron{trail,i}, C_coding_neuron{trail,i});
        overlap_BD(trail,i) = coding_overlap(B_coding_neuron{trail,i}, D_coding_neuron{trail,i});
        overlap_BE(trail,i) = coding_overlap(B_coding_neuron{trail,i}, E_coding_neuron{trail,i});
        overlap_BF(trail,i) = coding_overlap(B_coding_neuron{trail,i}, F_coding_neuron{trail,i});

        overlap_CD(trail,i) = coding_overlap(C_coding_neuron{trail,i}, D_coding_neuron{trail,i});
        overlap_CE(trail,i) = coding_overlap(C_coding_neuron{trail,i}, E_coding_neuron{trail,i});
        overlap_CF(trail,i) = coding_overlap(C_coding_neuron{trail,i}, F_coding_neuron{trail,i});

        overlap_DE(trail,i) = coding_overlap(D_coding_neuron{trail,i}, E_coding_neuron{trail,i});
        overlap_DF(trail,i) = coding_overlap(D_coding_neuron{trail,i}, F_coding_neuron{trail,i});

        overlap_EF(trail,i) = coding_overlap(E_coding_neuron{trail,i}, F_coding_neuron{trail,i});

        weight_A2A(trail,i) = mean_weight(exc_exc_weight_save,A_coding_neuron{trail,i},A_coding_neuron{trail,i},pos_S1(i));
        weight_A2B(trail,i) = mean_weight(exc_exc_weight_save,A_coding_neuron{trail,i},B_coding_neuron{trail,i},pos_S1(i));
        weight_A2C(trail,i) = mean_weight(exc_exc_weight_save,A_coding_neuron{trail,i},C_coding_neuron{trail,i},pos_S1(i));
        weight_A2D(trail,i) = mean_weight(exc_exc_weight_save,A_coding_neuron{trail,i},D_coding_neuron{trail,i},pos_S1(i));
        weight_A2E(trail,i) = mean_weight(exc_exc_weight_save,A_coding_neuron{trail,i},E_coding_neuron{trail,i},pos_S1(i));
        weight_A2F(trail,i) = mean_weight(exc_exc_weight_save,A_coding_neuron{trail,i},F_coding_neuron{trail,i},pos_S1(i));
        weight_A2O(trail,i) = mean_weight(exc_exc_weight_save,A_coding_neuron{trail,i},O_coding_neuron{trail,i},pos_S1(i));

        weight_B2A(trail,i) = mean_weight(exc_exc_weight_save,B_coding_neuron{trail,i},A_coding_neuron{trail,i},pos_S2(i));
        weight_B2B(trail,i) = mean_weight(exc_exc_weight_save,B_coding_neuron{trail,i},B_coding_neuron{trail,i},pos_S2(i));
        weight_B2C(trail,i) = mean_weight(exc_exc_weight_save,B_coding_neuron{trail,i},C_coding_neuron{trail,i},pos_S2(i));
        weight_B2D(trail,i) = mean_weight(exc_exc_weight_save,B_coding_neuron{trail,i},D_coding_neuron{trail,i},pos_S2(i));
        weight_B2E(trail,i) = mean_weight(exc_exc_weight_save,B_coding_neuron{trail,i},E_coding_neuron{trail,i},pos_S2(i));
        weight_B2F(trail,i) = mean_weight(exc_exc_weight_save,B_coding_neuron{trail,i},F_coding_neuron{trail,i},pos_S2(i));
        weight_B2O(trail,i) = mean_weight(exc_exc_weight_save,B_coding_neuron{trail,i},O_coding_neuron{trail,i},pos_S2(i));

        weight_C2A(trail,i) = mean_weight(exc_exc_weight_save,C_coding_neuron{trail,i},A_coding_neuron{trail,i},pos_S3(i));
        weight_C2B(trail,i) = mean_weight(exc_exc_weight_save,C_coding_neuron{trail,i},B_coding_neuron{trail,i},pos_S3(i));
        weight_C2C(trail,i) = mean_weight(exc_exc_weight_save,C_coding_neuron{trail,i},C_coding_neuron{trail,i},pos_S3(i));
        weight_C2D(trail,i) = mean_weight(exc_exc_weight_save,C_coding_neuron{trail,i},D_coding_neuron{trail,i},pos_S3(i));
        weight_C2E(trail,i) = mean_weight(exc_exc_weight_save,C_coding_neuron{trail,i},E_coding_neuron{trail,i},pos_S3(i));
        weight_C2F(trail,i) = mean_weight(exc_exc_weight_save,C_coding_neuron{trail,i},F_coding_neuron{trail,i},pos_S3(i));
        weight_C2O(trail,i) = mean_weight(exc_exc_weight_save,C_coding_neuron{trail,i},O_coding_neuron{trail,i},pos_S3(i));

        weight_D2A(trail,i) = mean_weight(exc_exc_weight_save,D_coding_neuron{trail,i},A_coding_neuron{trail,i},pos_S4(i));
        weight_D2B(trail,i) = mean_weight(exc_exc_weight_save,D_coding_neuron{trail,i},B_coding_neuron{trail,i},pos_S4(i));
        weight_D2C(trail,i) = mean_weight(exc_exc_weight_save,D_coding_neuron{trail,i},C_coding_neuron{trail,i},pos_S4(i));
        weight_D2D(trail,i) = mean_weight(exc_exc_weight_save,D_coding_neuron{trail,i},D_coding_neuron{trail,i},pos_S4(i));
        weight_D2E(trail,i) = mean_weight(exc_exc_weight_save,D_coding_neuron{trail,i},E_coding_neuron{trail,i},pos_S4(i));
        weight_D2F(trail,i) = mean_weight(exc_exc_weight_save,D_coding_neuron{trail,i},F_coding_neuron{trail,i},pos_S4(i));
        weight_D2O(trail,i) = mean_weight(exc_exc_weight_save,D_coding_neuron{trail,i},O_coding_neuron{trail,i},pos_S4(i));

        weight_E2A(trail,i) = mean_weight(exc_exc_weight_save,E_coding_neuron{trail,i},A_coding_neuron{trail,i},pos_S5(i));
        weight_E2B(trail,i) = mean_weight(exc_exc_weight_save,E_coding_neuron{trail,i},B_coding_neuron{trail,i},pos_S5(i));
        weight_E2C(trail,i) = mean_weight(exc_exc_weight_save,E_coding_neuron{trail,i},C_coding_neuron{trail,i},pos_S5(i));
        weight_E2D(trail,i) = mean_weight(exc_exc_weight_save,E_coding_neuron{trail,i},D_coding_neuron{trail,i},pos_S5(i));
        weight_E2E(trail,i) = mean_weight(exc_exc_weight_save,E_coding_neuron{trail,i},E_coding_neuron{trail,i},pos_S5(i));
        weight_E2F(trail,i) = mean_weight(exc_exc_weight_save,E_coding_neuron{trail,i},F_coding_neuron{trail,i},pos_S5(i));
        weight_E2O(trail,i) = mean_weight(exc_exc_weight_save,E_coding_neuron{trail,i},O_coding_neuron{trail,i},pos_S5(i));

        weight_F2A(trail,i) = mean_weight(exc_exc_weight_save,F_coding_neuron{trail,i},A_coding_neuron{trail,i},pos_S6(i));
        weight_F2B(trail,i) = mean_weight(exc_exc_weight_save,F_coding_neuron{trail,i},B_coding_neuron{trail,i},pos_S6(i));
        weight_F2C(trail,i) = mean_weight(exc_exc_weight_save,F_coding_neuron{trail,i},C_coding_neuron{trail,i},pos_S6(i));
        weight_F2D(trail,i) = mean_weight(exc_exc_weight_save,F_coding_neuron{trail,i},D_coding_neuron{trail,i},pos_S6(i));
        weight_F2E(trail,i) = mean_weight(exc_exc_weight_save,F_coding_neuron{trail,i},E_coding_neuron{trail,i},pos_S6(i));
        weight_F2F(trail,i) = mean_weight(exc_exc_weight_save,F_coding_neuron{trail,i},F_coding_neuron{trail,i},pos_S6(i));
        weight_F2O(trail,i) = mean_weight(exc_exc_weight_save,F_coding_neuron{trail,i},O_coding_neuron{trail,i},pos_S6(i));

        weight_O2A(trail,i) = mean_weight(exc_exc_weight_save,O_coding_neuron{trail,i},A_coding_neuron{trail,i},pos_S6(i));
        weight_O2B(trail,i) = mean_weight(exc_exc_weight_save,O_coding_neuron{trail,i},B_coding_neuron{trail,i},pos_S6(i));
        weight_O2C(trail,i) = mean_weight(exc_exc_weight_save,O_coding_neuron{trail,i},C_coding_neuron{trail,i},pos_S6(i));
        weight_O2D(trail,i) = mean_weight(exc_exc_weight_save,O_coding_neuron{trail,i},D_coding_neuron{trail,i},pos_S6(i));
        weight_O2E(trail,i) = mean_weight(exc_exc_weight_save,O_coding_neuron{trail,i},E_coding_neuron{trail,i},pos_S6(i));
        weight_O2F(trail,i) = mean_weight(exc_exc_weight_save,O_coding_neuron{trail,i},F_coding_neuron{trail,i},pos_S6(i));
        weight_O2O(trail,i) = mean_weight(exc_exc_weight_save,O_coding_neuron{trail,i},O_coding_neuron{trail,i},pos_S6(i));
    end
end

coding_neuron_num = (A_coding_neuron_num + B_coding_neuron_num + C_coding_neuron_num + D_coding_neuron_num + E_coding_neuron_num + F_coding_neuron_num)/single_pattern_num;
mean_coding_neuron_num = mean(coding_neuron_num);
std_coding_neuron_num = std(coding_neuron_num);

coding_neuron_fr = (A_coding_neuron_fr + B_coding_neuron_fr + C_coding_neuron_fr + D_coding_neuron_fr + E_coding_neuron_fr + F_coding_neuron_fr)/single_pattern_num;
mean_coding_neuron_fr = mean(coding_neuron_fr);
std_coding_neuron_fr = std(coding_neuron_fr);

coding_neuron{1,1} = A_coding_neuron;
coding_neuron{2,1} = B_coding_neuron;
coding_neuron{3,1} = C_coding_neuron;
coding_neuron{4,1} = D_coding_neuron;
coding_neuron{5,1} = E_coding_neuron;
coding_neuron{6,1} = F_coding_neuron;
coding_neuron{7,1} = O_coding_neuron;

save([Save_data_dir,'AI_inp=',num2str(inp_neuron_num),'_',num2str(single_pattern_num),'P_',num2str(single_pattern_show_time_new),'S_',...
    num2str(asso_pattern_show_time_new),'A_coding_neuron_inp-E=',num2str(inp_exc_ini),...
    '_E-E=',num2str(exc_exc_ini),'_E-I=',num2str(exc_inh_ini),'_I-E=',num2str(inh_exc_ini),'_I-I=',num2str(inh_inh_ini),'_pattern_time=',num2str(pattern_time),'_interval_time=',num2str(interval_time),...
    '_LTP_lr=',num2str(A2_plus_exc_novel),'_LTD_lr=',num2str(A2_minus_exc_novel),'_exc_ref=',num2str(exc_ref_mu),'_coding_thre=',num2str(coding_thre),'_trail_num=',num2str(trail_num),'.mat'],'coding_neuron');

save([Save_data_dir,'AI_inp=',num2str(inp_neuron_num),'_',num2str(single_pattern_num),'P_',num2str(single_pattern_show_time_new),'S_',...
    num2str(asso_pattern_show_time_new),'A_coding_neuron_num_inp-E=',num2str(inp_exc_ini),...
    '_E-E=',num2str(exc_exc_ini),'_E-I=',num2str(exc_inh_ini),'_I-E=',num2str(inh_exc_ini),'_I-I=',num2str(inh_inh_ini),'_pattern_time=',num2str(pattern_time),'_interval_time=',num2str(interval_time),...
    '_LTP_lr=',num2str(A2_plus_exc_novel),'_LTD_lr=',num2str(A2_minus_exc_novel),'_exc_ref=',num2str(exc_ref_mu),'_coding_thre=',num2str(coding_thre),'_trail_num=',num2str(trail_num),'.mat'],'coding_neuron_num');

save([Save_data_dir,'AI_inp=',num2str(inp_neuron_num),'_',num2str(single_pattern_num),'P_',num2str(single_pattern_show_time_new),'S_',...
    num2str(asso_pattern_show_time_new),'A_coding_neuron_fr_inp-E=',num2str(inp_exc_ini),...
    '_E-E=',num2str(exc_exc_ini),'_E-I=',num2str(exc_inh_ini),'_I-E=',num2str(inh_exc_ini),'_I-I=',num2str(inh_inh_ini),'_pattern_time=',num2str(pattern_time),'_interval_time=',num2str(interval_time),...
    '_LTP_lr=',num2str(A2_plus_exc_novel),'_LTD_lr=',num2str(A2_minus_exc_novel),'_exc_ref=',num2str(exc_ref_mu),'_coding_thre=',num2str(coding_thre),'_trail_num=',num2str(trail_num),'.mat'],'coding_neuron_fr');

overlap{1,1} = overlap_AB;
overlap{2,1} = overlap_AC;
overlap{3,1} = overlap_AD;
overlap{4,1} = overlap_AE;
overlap{5,1} = overlap_AF;

overlap{6,1} = overlap_BC;
overlap{7,1} = overlap_BD;
overlap{8,1} = overlap_BE;
overlap{9,1} = overlap_BF;

overlap{10,1} = overlap_CD;
overlap{11,1} = overlap_CE;
overlap{12,1} = overlap_CF;

overlap{13,1} = overlap_DE;
overlap{14,1} = overlap_DF;

overlap{15,1} = overlap_EF;

save([Save_data_dir,'AI_inp=',num2str(inp_neuron_num),'_',num2str(single_pattern_num),'P_',num2str(single_pattern_show_time_new),'S_',...
    num2str(asso_pattern_show_time_new),'A_overlap_inp-E=',num2str(inp_exc_ini),...
    '_E-E=',num2str(exc_exc_ini),'_E-I=',num2str(exc_inh_ini),'_I-E=',num2str(inh_exc_ini),'_I-I=',num2str(inh_inh_ini),'_pattern_time=',num2str(pattern_time),'_interval_time=',num2str(interval_time),...
    '_LTP_lr=',num2str(A2_plus_exc_novel),'_LTD_lr=',num2str(A2_minus_exc_novel),'_exc_ref=',num2str(exc_ref_mu),'_coding_thre=',num2str(coding_thre),'_trail_num=',num2str(trail_num),'.mat'],'overlap');

weight{1,1} = weight_A2A;
weight{2,1} = weight_A2B;
weight{3,1} = weight_A2C;
weight{4,1} = weight_A2D;
weight{5,1} = weight_A2E;
weight{6,1} = weight_A2F;
weight{7,1} = weight_A2O;

weight{1,2} = weight_B2A;
weight{2,2} = weight_B2B;
weight{3,2} = weight_B2C;
weight{4,2} = weight_B2D;
weight{5,2} = weight_B2E;
weight{6,2} = weight_B2F;
weight{7,2} = weight_B2O;

weight{1,3} = weight_C2A;
weight{2,3} = weight_C2B;
weight{3,3} = weight_C2C;
weight{4,3} = weight_C2D;
weight{5,3} = weight_C2E;
weight{6,3} = weight_C2F;
weight{7,3} = weight_C2O;

weight{1,4} = weight_D2A;
weight{2,4} = weight_D2B;
weight{3,4} = weight_D2C;
weight{4,4} = weight_D2D;
weight{5,4} = weight_D2E;
weight{6,4} = weight_D2F;
weight{7,4} = weight_D2O;

weight{1,5} = weight_E2A;
weight{2,5} = weight_E2B;
weight{3,5} = weight_E2C;
weight{4,5} = weight_E2D;
weight{5,5} = weight_E2E;
weight{6,5} = weight_E2F;
weight{7,5} = weight_E2O;

weight{1,6} = weight_F2A;
weight{2,6} = weight_F2B;
weight{3,6} = weight_F2C;
weight{4,6} = weight_F2D;
weight{5,6} = weight_F2E;
weight{6,6} = weight_F2F;
weight{7,6} = weight_F2O;

weight{1,7} = weight_O2A;
weight{2,7} = weight_O2B;
weight{3,7} = weight_O2C;
weight{4,7} = weight_O2D;
weight{5,7} = weight_O2E;
weight{6,7} = weight_O2F;
weight{7,7} = weight_O2O;

save([Save_data_dir,'AI_inp=',num2str(inp_neuron_num),'_',num2str(single_pattern_num),'P_',num2str(single_pattern_show_time_new),'S_',...
    num2str(asso_pattern_show_time_new),'A_weight_inp-E=',num2str(inp_exc_ini),...
    '_E-E=',num2str(exc_exc_ini),'_E-I=',num2str(exc_inh_ini),'_I-E=',num2str(inh_exc_ini),'_I-I=',num2str(inh_inh_ini),'_pattern_time=',num2str(pattern_time),'_interval_time=',num2str(interval_time),...
    '_LTP_lr=',num2str(A2_plus_exc_novel),'_LTD_lr=',num2str(A2_minus_exc_novel),'_exc_ref=',num2str(exc_ref_mu),'_coding_thre=',num2str(coding_thre),'_trail_num=',num2str(trail_num),'.mat'],'weight');