%%% Plot neural coding overlap and coding neuron number across time
%%% Author: Ye Wang

close all
clear
clc

today_date = datestr(now,'yymmdd');

main_dir = 'D:\Memory_integration\';

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
    num2str(asso_pattern_show_time_new),'A_overlap_inp-E=',num2str(inp_exc_ini),...
    '_E-E=',num2str(exc_exc_ini),'_E-I=',num2str(exc_inh_ini),'_I-E=',num2str(inh_exc_ini),'_I-I=',num2str(inh_inh_ini),'_pattern_time=',num2str(pattern_time),'_interval_time=',num2str(interval_time),...
    '_LTP_lr=',num2str(A2_plus_exc_novel),'_LTD_lr=',num2str(A2_minus_exc_novel),'_exc_ref=',num2str(exc_ref_mu),'_coding_thre=',num2str(coding_thre),'_trail_num=',num2str(trail_num),'.mat']);

load([Save_data_dir,'AI_inp=',num2str(inp_neuron_num),'_',num2str(single_pattern_num),'P_',num2str(single_pattern_show_time_new),'S_',...
    num2str(asso_pattern_show_time_new),'A_coding_neuron_num_inp-E=',num2str(inp_exc_ini),...
    '_E-E=',num2str(exc_exc_ini),'_E-I=',num2str(exc_inh_ini),'_I-E=',num2str(inh_exc_ini),'_I-I=',num2str(inh_inh_ini),'_pattern_time=',num2str(pattern_time),'_interval_time=',num2str(interval_time),...
    '_LTP_lr=',num2str(A2_plus_exc_novel),'_LTD_lr=',num2str(A2_minus_exc_novel),'_exc_ref=',num2str(exc_ref_mu),'_coding_thre=',num2str(coding_thre),'_trail_num=',num2str(trail_num),'.mat']);

load([Save_data_dir,'AI_inp=',num2str(inp_neuron_num),'_',num2str(single_pattern_num),'P_',num2str(single_pattern_show_time_new),'S_',...
    num2str(asso_pattern_show_time_new),'A_total_coding_neuron_num_inp-E=',num2str(inp_exc_ini),...
    '_E-E=',num2str(exc_exc_ini),'_E-I=',num2str(exc_inh_ini),'_I-E=',num2str(inh_exc_ini),'_I-I=',num2str(inh_inh_ini),'_pattern_time=',num2str(pattern_time),'_interval_time=',num2str(interval_time),...
    '_LTP_lr=',num2str(A2_plus_exc_novel),'_LTD_lr=',num2str(A2_minus_exc_novel),'_exc_ref=',num2str(exc_ref_mu),'_coding_thre=',num2str(coding_thre),'_trail_num=',num2str(trail_num),'.mat']);


%%% Overlap
overlap_AB = overlap{1,1};
overlap_AC = overlap{2,1};
overlap_AD = overlap{3,1};
overlap_AE = overlap{4,1};
overlap_AF = overlap{5,1};

overlap_BC = overlap{6,1};
overlap_BD = overlap{7,1};
overlap_BE = overlap{8,1};
overlap_BF = overlap{9,1};

overlap_CD = overlap{10,1};
overlap_CE = overlap{11,1};
overlap_CF = overlap{12,1};

overlap_DE = overlap{13,1};
overlap_DF = overlap{14,1};

overlap_EF = overlap{15,1};

overlap_lag_0 = (overlap_AB + overlap_BC + overlap_DE + overlap_EF)/4;
overlap_lag_1 = (overlap_AC + overlap_DF)/2;
overlap_lag_inf = (overlap_AD + overlap_AE + overlap_AF + overlap_BD + overlap_BE + overlap_BF + overlap_CD + overlap_CE + overlap_CF)/9;

xlswrite([Save_data_dir,'overlap_lag_0.xlsx'], overlap_lag_0);
xlswrite([Save_data_dir,'overlap_lag_1.xlsx'], overlap_lag_1);
xlswrite([Save_data_dir,'overlap_lag_inf.xlsx'], overlap_lag_inf);


mean_overlap_lag_0 = mean(overlap_lag_0);
mean_overlap_lag_1 = mean(overlap_lag_1);
mean_overlap_lag_inf = mean(overlap_lag_inf);

std_overlap_lag_0 = std(overlap_lag_0);
std_overlap_lag_1 = std(overlap_lag_1);
std_overlap_lag_inf = std(overlap_lag_inf);

X = [1:total_show_time,total_show_time:-1:1];

fill_lag_0_y1 = mean_overlap_lag_0 - std_overlap_lag_0;
fill_lag_0_y2 = mean_overlap_lag_0 + std_overlap_lag_0;
fill_lag_0_Y = [fill_lag_0_y1,fliplr(fill_lag_0_y2)];

fill_lag_1_y1 = mean_overlap_lag_1 - std_overlap_lag_1;
fill_lag_1_y2 = mean_overlap_lag_1 + std_overlap_lag_1;
fill_lag_1_Y = [fill_lag_1_y1,fliplr(fill_lag_1_y2)];

fill_lag_inf_y1 = mean_overlap_lag_inf - std_overlap_lag_inf;
fill_lag_inf_y2 = mean_overlap_lag_inf + std_overlap_lag_inf;
fill_lag_inf_Y = [fill_lag_inf_y1,fliplr(fill_lag_inf_y2)];

color_list = [0.85,0.33,0.1;
    0.00,0.45,0.74;
    0.93,0.69,0.13];
figure,
set(gcf,'Position',[129.4000 255.4000 1128 420])
h1 = fill(X,fill_lag_0_Y,color_list(1,:));
set(h1,'edgealpha',0,'facealpha',0.1) 
hold on
p1 = errorbar(1:total_show_time,mean_overlap_lag_0,std_overlap_lag_0,':','LineWidth',2,'Color',color_list(1,:));
hold on
h2 = fill(X,fill_lag_1_Y,color_list(2,:));
set(h2,'edgealpha',0,'facealpha',0.1) 
hold on
p2 = errorbar(1:total_show_time,mean_overlap_lag_1,std_overlap_lag_1,':','LineWidth',2,'Color',color_list(2,:));
hold on
h3 = fill(X,fill_lag_inf_Y,color_list(3,:));
set(h3,'edgealpha',0,'facealpha',0.1) 
p3 = errorbar(1:total_show_time,mean_overlap_lag_inf,std_overlap_lag_inf,':','LineWidth',2,'Color',color_list(3,:));
hold off
set(gca,'XTick',1:total_show_time);
set(gca,'XTickLabel',{'T1','T2','T3','T4','T5','T6','T7','T8','T9','T10','T11','T12','T13','T14','T15','T16','T17','T18','T19','T20'}); 
legend([p1 p2 p3],{'0бу','1бу','inf'},'Location','NorthWest')
xlabel('Single memory item show time')
ylabel('Neural coding overlap')
xlim([0.5,total_show_time+0.5])
box off

%%% Coding neuron number
mean_coding_neuron_num = mean(coding_neuron_num);
std_coding_neuron_num = std(coding_neuron_num);
xlswrite([Save_data_dir,'coding_neuron_num.xlsx'], coding_neuron_num);

mean_total_coding_neuron_num = mean(total_coding_neuron_num);
std_total_coding_neuron_num = std(total_coding_neuron_num);

figure,
set(gcf,'Position',[211.8000 182.2000 1.1264e+03 289.6000])
b1 = bar(mean_total_coding_neuron_num,0.5);
hold on
errorbar(1:total_show_time,mean_total_coding_neuron_num,std_total_coding_neuron_num, 'k','Linestyle', 'None');
hold on
b2 = bar(mean_coding_neuron_num,0.5);
hold on
errorbar(1:total_show_time,mean_coding_neuron_num,std_coding_neuron_num, 'k','Linestyle', 'None');

legend([b1 b2],{'total','single'},'Location','NorthEast')
set(gca,'XTick',1:total_show_time);
set(gca,'XTickLabel',{'T1','T2','T3','T4','T5','T6','T7','T8','T9','T10','T11','T12','T13','T14','T15','T16','T17','T18','T19','T20'}); 
xlabel('Single memory item show time')
ylabel('Coding neuron number')
xlim([0.5,total_show_time+0.5])
box off