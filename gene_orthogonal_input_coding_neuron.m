%%% Generate othorgonal coding neuron for input layer
%%% The pattern number can be defined by user
%%% The coding neuron for each trail is different
%%% Author: Ye Wang

close all
clear 
clc

if ~exist('D:\Memory_integration\','dir')==1
   mkdir('D:\Memory_integration\');
end

main_dir = 'D:\Memory_integration\';

if ~exist([main_dir,'stim\'],'dir')==1
   mkdir([main_dir,'stim\']);
end

inp_neuron_num = 400; 
sparsity = 0.1; 
coding_neuron_num = sparsity * inp_neuron_num;
pattern_num = 6;
trail_num = 30;


for trail = 1:trail_num
    shuffle_index = randperm(inp_neuron_num);
    for i = 1:pattern_num
        pattern(i,:) = shuffle_index(coding_neuron_num*(i-1)+1:coding_neuron_num*i);
    end
    save([main_dir,'stim\','shuffle_inp=',num2str(inp_neuron_num),'_pattern_num=',num2str(pattern_num),'_trail=',num2str(trail),'.mat'],'pattern');
end

    




