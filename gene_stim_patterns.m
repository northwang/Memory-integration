%%% Read the input coding neuron file and generate stimuli file
%%% Author: Ye Wang

close all
clear 
clc

main_dir = 'D:\Memory_integration\stim\';
show_pattern_num = 6;
single_pattern_show_time = 10; 
asso_group_num = 4;  
show_asso_group_num = 4;
asso_group_show_time = 10; 
validation_show_time = 10; 
pattern_time = 1000;
interval_time = 1000;

total_stim_length = show_pattern_num*single_pattern_show_time + show_asso_group_num*asso_group_show_time + show_pattern_num*validation_show_time;

inp_neuron_num = 400;

trail_num = 30;

for trail = 1:trail_num
    
    shuffle_list = 1:total_stim_length-1;

    for k = 1:single_pattern_show_time 
        shuffle_list((k-1)*show_pattern_num + 1:k*show_pattern_num) = (k-1)*show_pattern_num + randperm(show_pattern_num);
    end

    for k = 1:asso_group_show_time
        shuffle_list(show_pattern_num*single_pattern_show_time + (k-1)*(show_pattern_num + asso_group_num) + 1:show_pattern_num*single_pattern_show_time + (k-1)*(show_pattern_num + asso_group_num) + asso_group_num) =...
            show_pattern_num*single_pattern_show_time + (k-1)*(show_pattern_num + asso_group_num) + randperm(asso_group_num);
    end
    
    for k = 1:asso_group_show_time
        shuffle_list(show_pattern_num*single_pattern_show_time + (k-1)*(show_pattern_num + asso_group_num) + asso_group_num + 1:show_pattern_num*single_pattern_show_time + k*(show_pattern_num + asso_group_num)) =...
            show_pattern_num*single_pattern_show_time + (k-1)*(show_pattern_num + asso_group_num) + asso_group_num + randperm(show_pattern_num);
    end
    
    shuffle_list = shuffle_list';

    
    coding_file = ['shuffle_inp=',num2str(inp_neuron_num),'_pattern_num=',num2str(show_pattern_num),'_trail=',num2str(trail)];

    coding_path = [main_dir,coding_file,'.mat'];
    load(coding_path);

    % Associated pattern
    asso_pattern{1,:} = [pattern(1,:),pattern(2,:)]; 
    asso_pattern{2,:} = [pattern(2,:),pattern(3,:)]; 
    asso_pattern{3,:} = [pattern(4,:),pattern(5,:)];
    asso_pattern{4,:} = [pattern(5,:),pattern(6,:)]; 
    
    stim_list = cell(1,3);
    for i = 1:total_stim_length
        stim_list{2*i-1,1} = 'blank';
        stim_list{2*i-1,2} = interval_time;
        stim_list{2*i-1,3} = 0;
        if i <= show_pattern_num*single_pattern_show_time 
            if rem(shuffle_list(i),show_pattern_num) ~= 0 
                stim_list{2*i,1} = ['S',num2str(rem(shuffle_list(i),show_pattern_num))];
                stim_list{2*i,3} = pattern(rem(shuffle_list(i),show_pattern_num),:);
            else
                stim_list{2*i,1} = ['S',num2str(show_pattern_num)];
                stim_list{2*i,3} = pattern(show_pattern_num,:);
            end
        else
            if rem(shuffle_list(i)-show_pattern_num*single_pattern_show_time,asso_group_num+show_pattern_num) == 0
                stim_list{2*i,1} = ['S',num2str(show_pattern_num)];
                stim_list{2*i,3} = pattern(show_pattern_num,:);
            elseif (rem(i-show_pattern_num*single_pattern_show_time,asso_group_num+show_pattern_num) >= 1) &&...
                    (rem(i-show_pattern_num*single_pattern_show_time,asso_group_num+show_pattern_num) <= asso_group_num) 
                stim_list{2*i,1} = ['A',num2str(rem(shuffle_list(i)-show_pattern_num*single_pattern_show_time,asso_group_num+show_pattern_num))];
                stim_list{2*i,3} = asso_pattern{rem(shuffle_list(i)-show_pattern_num*single_pattern_show_time,asso_group_num+show_pattern_num),:};
            else 
                stim_list{2*i,1} = ['S',num2str(rem(shuffle_list(i)-show_pattern_num*single_pattern_show_time,asso_group_num+show_pattern_num) - asso_group_num)];
                stim_list{2*i,3} = pattern(rem(shuffle_list(i)-show_pattern_num*single_pattern_show_time,asso_group_num+show_pattern_num) - asso_group_num,:);
            end
        end
        stim_list{2*i,2} = pattern_time;
    end

    save([main_dir,'AI_shuffle_show_inp=',num2str(inp_neuron_num),'_pattern_num=',num2str(show_pattern_num),'_single_reps=',num2str(single_pattern_show_time),'_asso_reps=',num2str(asso_group_show_time)...
        '_test_reps=',num2str(validation_show_time),'_pattern_time=',num2str(pattern_time),'_interval_time=',num2str(interval_time),'_trail=',num2str(trail),'.mat'],'stim_list'); 
end


