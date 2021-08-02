%%% Plot weight change of pattern pairs with different degree
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
exc_exc_conn_rate = 0.2;
single_pattern_num = 6;
single_pattern_show_time_new = 10; 
single_pattern_show_time_old = 0;
asso_pattern_show_time_new = 10;
asso_pattern_show_time_old = 0;
validation_show_time = asso_pattern_show_time_new;
pattern_show_time = single_pattern_show_time_new; 
total_show_time = single_pattern_show_time_new + validation_show_time;
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

trail_num = 30;

load([Save_data_dir,'AI_inp=',num2str(inp_neuron_num),'_',num2str(single_pattern_num),'P_',num2str(single_pattern_show_time_new),'S_',...
    num2str(asso_pattern_show_time_new),'A_weight_inp-E=',num2str(inp_exc_ini),...
    '_E-E=',num2str(exc_exc_ini),'_E-I=',num2str(exc_inh_ini),'_I-E=',num2str(inh_exc_ini),'_I-I=',num2str(inh_inh_ini),'_pattern_time=',num2str(pattern_time),'_interval_time=',num2str(interval_time),...
    '_LTP_lr=',num2str(A2_plus_exc_novel),'_LTD_lr=',num2str(A2_minus_exc_novel),'_exc_ref=',num2str(exc_ref_mu),'_coding_thre=',num2str(coding_thre),'_trail_num=',num2str(trail_num),'.mat']);

load([Save_data_dir,'AI_inp=',num2str(inp_neuron_num),'_',num2str(single_pattern_num),'P_',num2str(single_pattern_show_time_new),'S_',...
    num2str(asso_pattern_show_time_new),'A_coding_neuron_fr_inp-E=',num2str(inp_exc_ini),...
    '_E-E=',num2str(exc_exc_ini),'_E-I=',num2str(exc_inh_ini),'_I-E=',num2str(inh_exc_ini),'_I-I=',num2str(inh_inh_ini),'_pattern_time=',num2str(pattern_time),'_interval_time=',num2str(interval_time),...
    '_LTP_lr=',num2str(A2_plus_exc_novel),'_LTD_lr=',num2str(A2_minus_exc_novel),'_exc_ref=',num2str(exc_ref_mu),'_coding_thre=',num2str(coding_thre),'_trail_num=',num2str(trail_num),'.mat']);

mean_coding_neuron_fr = mean(coding_neuron_fr);
std_coding_neuron_fr = std(coding_neuron_fr);

figure,
set(gcf,'Position',[211.8000 182.2000 1.1264e+03 289.6000])
b1 = bar(mean_coding_neuron_fr,0.5);
hold on
errorbar(1:total_show_time,mean_coding_neuron_fr,std_coding_neuron_fr, 'k','Linestyle', 'None');
hold off
set(gca,'XTick',1:total_show_time);
set(gca,'XTickLabel',{'T1','T2','T3','T4','T5','T6','T7','T8','T9','T10','T11','T12','T13','T14','T15','T16','T17','T18','T19','T20'}); 
xlabel('Single memory item show time')
ylabel('Coding neuron firing rate (spks/s)')
xlim([0.5,total_show_time+0.5])
box off

weight_A2A = weight{1,1};
weight_A2B = weight{2,1};
weight_A2C = weight{3,1};
weight_A2D = weight{4,1};
weight_A2E = weight{5,1};
weight_A2F = weight{6,1};
weight_A2O = weight{7,1};

weight_B2A = weight{1,2};
weight_B2B = weight{2,2};
weight_B2C = weight{3,2};
weight_B2D = weight{4,2};
weight_B2E = weight{5,2};
weight_B2F = weight{6,2};
weight_B2O = weight{7,2};

weight_C2A = weight{1,3};
weight_C2B = weight{2,3};
weight_C2C = weight{3,3};
weight_C2D = weight{4,3};
weight_C2E = weight{5,3};
weight_C2F = weight{6,3};
weight_C2O = weight{7,3};

weight_D2A = weight{1,4};
weight_D2B = weight{2,4};
weight_D2C = weight{3,4};
weight_D2D = weight{4,4};
weight_D2E = weight{5,4};
weight_D2F = weight{6,4};
weight_D2O = weight{7,4};

weight_E2A = weight{1,5};
weight_E2B = weight{2,5};
weight_E2C = weight{3,5};
weight_E2D = weight{4,5};
weight_E2E = weight{5,5};
weight_E2F = weight{6,5};
weight_E2O = weight{7,5};

weight_F2A = weight{1,6};
weight_F2B = weight{2,6};
weight_F2C = weight{3,6};
weight_F2D = weight{4,6};
weight_F2E = weight{5,6};
weight_F2F = weight{6,6};
weight_F2O = weight{7,6};

weight_O2A = weight{1,7};
weight_O2B = weight{2,7};
weight_O2C = weight{3,7};
weight_O2D = weight{4,7};
weight_O2E = weight{5,7};
weight_O2F = weight{6,7};
weight_O2O = weight{7,7};

weight_inner = (weight_A2A + weight_B2B + weight_C2C + weight_D2D + weight_E2E + weight_F2F)/6;
weight_lag_0 = (weight_A2B + weight_B2A + weight_B2C + weight_C2B + weight_D2E + weight_E2D + weight_E2F + weight_F2E)/8;
weight_lag_1 = (weight_A2C + weight_C2A + weight_D2F + weight_F2D)/4;
weight_lag_inf = (weight_A2D + weight_A2E + weight_A2F + ...
    weight_B2D + weight_B2E + weight_B2F + ...
    weight_C2D + weight_C2E + weight_C2F + ...
    weight_D2A + weight_D2B + weight_D2C + ...
    weight_E2A + weight_E2B + weight_E2C + ...
    weight_F2A + weight_F2B + weight_F2C)/18;
weight_other = (weight_A2O + weight_B2O + weight_C2O + weight_D2O + weight_E2O + weight_F2O + ...
    weight_O2A + weight_O2B + weight_O2C + weight_O2D + weight_O2E + weight_O2F)/12;

save([Save_data_dir,'AI_inp=',num2str(inp_neuron_num),'_',num2str(single_pattern_num),'P_',num2str(single_pattern_show_time_new),'S_',...
    num2str(asso_pattern_show_time_new),'A_weight_inp-E=',num2str(inp_exc_ini),...
    '_E-E=',num2str(exc_exc_ini),'_E-I=',num2str(exc_inh_ini),'_I-E=',num2str(inh_exc_ini),'_I-I=',num2str(inh_inh_ini),'_pattern_time=',num2str(pattern_time),'_interval_time=',num2str(interval_time),...
    '_LTP_lr=',num2str(A2_plus_exc_novel),'_LTD_lr=',num2str(A2_minus_exc_novel),'_exc_ref=',num2str(exc_ref_mu),'_coding_thre=',num2str(coding_thre),'_trail_num=',num2str(trail_num),'.mat'],'weight');

xlswrite([Save_data_dir,'weight_inner.xlsx'], weight_inner);
xlswrite([Save_data_dir,'weight_lag_0.xlsx'], weight_lag_0);
xlswrite([Save_data_dir,'weight_lag_1.xlsx'], weight_lag_1);
xlswrite([Save_data_dir,'weight_lag_inf.xlsx'], weight_lag_inf);
xlswrite([Save_data_dir,'weight_other.xlsx'], weight_other);

mean_weight_inner = mean(weight_inner);
mean_weight_lag_0 = mean(weight_lag_0);
mean_weight_lag_1 = mean(weight_lag_1);
mean_weight_lag_inf = mean(weight_lag_inf);
mean_weight_other = mean(weight_other);

std_weight_inner = std(weight_inner);
std_weight_lag_0 = std(weight_lag_0);
std_weight_lag_1 = std(weight_lag_1);
std_weight_lag_inf = std(weight_lag_inf);
std_weight_other = std(weight_other);

X = [1:total_show_time,total_show_time:-1:1];

fill_weight_inner_y1 = mean_weight_inner - std_weight_inner;
fill_weight_inner_y2 = mean_weight_inner + std_weight_inner;
fill_weight_inner_Y = [fill_weight_inner_y1,fliplr(fill_weight_inner_y2)];

fill_weight_lag_0_y1 = mean_weight_lag_0 - std_weight_lag_0;
fill_weight_lag_0_y2 = mean_weight_lag_0 + std_weight_lag_0;
fill_weight_lag_0_Y = [fill_weight_lag_0_y1,fliplr(fill_weight_lag_0_y2)];

fill_weight_lag_1_y1 = mean_weight_lag_1 - std_weight_lag_1;
fill_weight_lag_1_y2 = mean_weight_lag_1 + std_weight_lag_1;
fill_weight_lag_1_Y = [fill_weight_lag_1_y1,fliplr(fill_weight_lag_1_y2)];

fill_weight_lag_inf_y1 = mean_weight_lag_inf - std_weight_lag_inf;
fill_weight_lag_inf_y2 = mean_weight_lag_inf + std_weight_lag_inf;
fill_weight_lag_inf_Y = [fill_weight_lag_inf_y1,fliplr(fill_weight_lag_inf_y2)];

fill_weight_other_y1 = mean_weight_other - std_weight_other;
fill_weight_other_y2 = mean_weight_other + std_weight_other;
fill_weight_other_Y = [fill_weight_other_y1,fliplr(fill_weight_other_y2)];


color_list = [0.85,0.33,0.1;
    0.00,0.45,0.74;
    0.93,0.69,0.13;
    0.49,0.18,0.56;
    0.47,0.67,0.19];

figure,
set(gcf,'Position',[129.4000 255.4000 1128 420])
h1 = fill(X,fill_weight_inner_Y,color_list(1,:));
set(h1,'edgealpha',0,'facealpha',0.1) 
hold on
p1 = errorbar(1:total_show_time,mean_weight_inner,std_weight_inner,':','LineWidth',2,'Color',color_list(1,:));
hold on
h2 = fill(X,fill_weight_lag_0_Y,color_list(2,:));
set(h2,'edgealpha',0,'facealpha',0.1) 
hold on
p2 = errorbar(1:total_show_time,mean_weight_lag_0,std_weight_lag_0,':','LineWidth',2,'Color',color_list(2,:));
hold on
h3 = fill(X,fill_weight_lag_1_Y,color_list(3,:));
set(h3,'edgealpha',0,'facealpha',0.1) 
hold on
p3 = errorbar(1:total_show_time,mean_weight_lag_1,std_weight_lag_1,':','LineWidth',2,'Color',color_list(3,:));
hold on
h4 = fill(X,fill_weight_lag_inf_Y,color_list(4,:));
set(h4,'edgealpha',0,'facealpha',0.1) 
p4 = errorbar(1:total_show_time,mean_weight_lag_inf,std_weight_lag_inf,':','LineWidth',2,'Color',color_list(4,:));
hold on
h5 = fill(X,fill_weight_other_Y,color_list(5,:));
set(h5,'edgealpha',0,'facealpha',0.1) 
p5 = errorbar(1:total_show_time,mean_weight_other,std_weight_other,':','LineWidth',2,'Color',color_list(5,:));
hold off
set(gca,'XTick',1:total_show_time);
set(gca,'XTickLabel',{'T1','T2','T3','T4','T5','T6','T7','T8','T9','T10','T11','T12','T13','T14','T15','T16','T17','T18','T19','T20'}); 
legend([p1 p2 p3 p4 p5],{'inner','0бу','1бу','inf','others'},'Location','NorthEast')
xlabel('Single memory item show time')
ylabel('Weight')
xlim([0.5,total_show_time+0.5])
box off

