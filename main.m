%%% Formation of association with assemblies formed first
%%% Simulation with E->E triplet STDP and homeostatic plasticity
%%% Author: Ye Wang

close all
clear
clc

main_dir = 'D:\Memory_integration\';

if ~exist([main_dir,'figures\'],'dir')==1
   mkdir([main_dir,'figures\']);
end

if ~exist([main_dir,'data\'],'dir')==1
   mkdir([main_dir,'data\']);
end

Save_figure_dir = [main_dir,'figures\'];
Save_data_dir = [main_dir,'data\'];


%% Neuron numbers
inp_neuron_num = 400;
exc_neuron_num = 800; 
inh_neuron_num = 200; 

%% Load stimuli file
single_pattern_num = 6;
asso_pattern_num = 4;
single_pattern_show_time = 10;
asso_pattern_show_time = 10;
validation_show_time = 10;

pattern_time = 1000;
interval_time = 1000;
Nskip = 1000; % how often (in number of timesteps) to save weight

trail_num = 30;

for trail = 1:trail_num

    start_time = datetime('now');

    stim_file = ['stim\AI_shuffle_show_inp=',num2str(inp_neuron_num),'_pattern_num=',num2str(single_pattern_num),'_single_reps=',num2str(single_pattern_show_time),'_asso_reps=',num2str(asso_pattern_show_time)...
        '_test_reps=',num2str(validation_show_time),'_pattern_time=',num2str(pattern_time),'_interval_time=',num2str(interval_time),'_trail=',num2str(trail)]; 

    stim_path = [main_dir,stim_file,'.mat'];
    load(stim_path); 

    stim_time = sum(cell2mat(stim_list(:,2))) + 10; 

    %% Connect rate
    inp_exc_conn_rate = 0.2; 
    exc_exc_conn_rate = 0.2;
    exc_inh_conn_rate = 0.2;
    inh_exc_conn_rate = 0.2; 
    inh_inh_conn_rate = 0.2;
    inp_inh_conn_rate = 0.1;

    inp_exc_ini = 30; 
    exc_exc_ini = 10; 
    exc_inh_ini = 5;
    inh_exc_ini = 30; 
    inh_inh_ini = 5;
    inp_inh_ini = 20;

    max_exc_exc_weight = 10*exc_exc_ini;
    min_exc_exc_weight = 0.1*exc_exc_ini;


    %% Input parameters
    coding_fr = 200; % Input coding neuron firing rate
    background_fr = 0.1; 
    noise_fr = 0.1; 

    % % Response kernel
    tau_f = 20;
    u_thre = 30;  % firing threshold

    % refractory period
    exc_ref_mu = 20;
    ref_shape = 2;

    inh_ref_mu = 4;

    exc_ref_max = 10*exc_ref_mu;
    inh_ref_max = 10*inh_ref_mu;

    %% Parameters for plasticity

    A2_plus_exc_novel = 0.4;
    A2_minus_exc_novel = 0.02;  
    A3_plus_exc_novel = 0.4;  
    A3_minus_exc_novel = 0.02; 

    tau_r1 = 25;
    tau_r2 = 25;
    tau_o1 = 1000;
    tau_o2 = 25;

    r1_exc_exc_0 = 0;
    r2_exc_exc_0 = 0;
    o1_exc_exc_0 = 0;
    o2_exc_exc_0 = 0;

    r1_exc_exc = zeros(exc_neuron_num,exc_neuron_num,2);
    r2_exc_exc = zeros(exc_neuron_num,exc_neuron_num,2);
    o1_exc_exc = zeros(exc_neuron_num,exc_neuron_num,2);
    o2_exc_exc = zeros(exc_neuron_num,exc_neuron_num,2);

    exc_exc_length = 6;

    inh_exc_length = 3;
    inh_exc_update_flag = zeros(inh_neuron_num,exc_neuron_num);

    %% parameters for kernel
    inp_res_ker_length = 6;
    exc_res_ker_length = 6;
    inh_res_ker_length = 3;

    %% input neurons
    inp_res_ker = zeros(inp_neuron_num,inp_res_ker_length);
    inp_spike_time = zeros(inp_neuron_num,1); 
    inp_spike_count = zeros(inp_neuron_num,1);
    inp_spike_flag = zeros(inp_neuron_num,1); 
    inp_last_spike = zeros(inp_neuron_num,1); 
    inp_fr = zeros(inp_neuron_num,1); 

    %% Exc and inh neurons
    exc_spike_time = zeros(exc_neuron_num,1); 
    exc_spike_count = zeros(exc_neuron_num,1); 
    exc_spike_flag = zeros(exc_neuron_num,1); 
    exc_last_spike = zeros(exc_neuron_num,1);
    exc_fr = zeros(exc_neuron_num,1); 

    exc_u = zeros(exc_neuron_num,1); 
    exc = zeros(exc_neuron_num,exc_ref_max); 
    exc_res_ker = zeros(exc_neuron_num,exc_res_ker_length);

    inh_spike_time = zeros(inh_neuron_num,1); 
    inh_spike_count = zeros(inh_neuron_num,1); 
    inh_spike_flag = zeros(inh_neuron_num,1); 
    inh_last_spike = zeros(inh_neuron_num,1); 
    inh_fr = zeros(inh_neuron_num,1); 

    inh_u = zeros(inh_neuron_num,1); 
    inh = zeros(inh_neuron_num,inh_ref_max);
    inh_res_ker = zeros(inh_neuron_num,inh_res_ker_length); 

    %% load input->E, E->E, E->I, I->E, I->I connections
    [inp_exc_weight,exc_exc_weight,exc_inh_weight,inh_exc_weight,inh_inh_weight,inp_inh_weight] = ...
        conn_weight_new(inp_neuron_num,exc_neuron_num,inh_neuron_num,inp_exc_conn_rate,exc_exc_conn_rate,exc_inh_conn_rate,inh_exc_conn_rate,inh_inh_conn_rate,inp_inh_conn_rate,...
        inp_exc_ini,exc_exc_ini,exc_inh_ini, inh_exc_ini, inh_inh_ini,inp_inh_ini);

    exc_exc_weight_save = reshape(repmat(exc_exc_weight,1,floor(stim_time/Nskip) + 1),exc_neuron_num,exc_neuron_num,floor(stim_time/Nskip) + 1); 
    exc_exc_update_flag = zeros(exc_neuron_num,exc_neuron_num);    
    exc_exc_updated_w = reshape(repmat(exc_exc_weight,1,exc_exc_length),exc_neuron_num,exc_neuron_num,exc_exc_length);
    exc_exc_weight_0 = sum(exc_exc_weight);
    exc_exc_weight_pre_num = exc_exc_weight_0/exc_exc_ini;
    exc_exc_weight_mask = exc_exc_weight/exc_exc_ini;

    fr_thre_noise = noise_fr*ones(inp_neuron_num,1);

    blank_show_time = 0;
    S_show_time = 0;
    A_show_time = 0;
    plasticity = 'OFF';

    t = 10;
    for stim = 1:size(stim_list,1)
        current_stim = char(stim_list(stim,1));
        if strcmp(current_stim,'blank')
            blank_show_time = blank_show_time + 1;
        elseif strcmp(current_stim(1),'S')
            S_show_time = S_show_time + 1/single_pattern_num;
            if A_show_time == 0
                plasticity = 'ON';
            else 
                plasticity = 'OFF';
            end
        elseif strcmp(current_stim(1),'A')
            A_show_time = A_show_time + 1/asso_pattern_num;
            plasticity = 'ON';
        end

        if strcmp(plasticity,'ON')
            A2_plus_exc = A2_plus_exc_novel; 
            A2_minus_exc = A2_minus_exc_novel; 
            A3_plus_exc = A3_plus_exc_novel; 
            A3_minus_exc = A3_minus_exc_novel; 
        else
            A2_plus_exc = 0; 
            A2_minus_exc = 0; 
            A3_plus_exc = 0; 
            A3_minus_exc = 0; 
        end 

        for step = 1:stim_list{stim,2}
            if rem(t,Nskip) == 0
                print_info = '%.4f completed，%s left\n';
                fprintf(print_info,t/stim_time,(datetime('now') - start_time)*(stim_time-t)/t) % print run time
                if t == Nskip
                    inp_fr(:,t/Nskip) = inp_spike_count*1000/Nskip;
                    exc_fr(:,t/Nskip) = exc_spike_count*1000/Nskip;
                    inh_fr(:,t/Nskip) = inh_spike_count*1000/Nskip;
                else
                    inp_fr(:,t/Nskip) = inp_spike_count*1000/Nskip - sum(inp_fr(:,1:(t/Nskip) - 1),2);
                    exc_fr(:,t/Nskip) = exc_spike_count*1000/Nskip - sum(exc_fr(:,1:(t/Nskip) - 1),2);
                    inh_fr(:,t/Nskip) = inh_spike_count*1000/Nskip - sum(inh_fr(:,1:(t/Nskip) - 1),2);
                end
            end

            t = t + 1;
            % input neurons
            inp = rand(inp_neuron_num,1); 

            if strcmp(stim_list(stim,1),'blank')
                inp(inp <= fr_thre_noise/1000) = 1;
            else
                fr_thre = background_fr*ones(inp_neuron_num,1);
                stim_pattern = stim_list{stim,3};
                fr_thre(stim_pattern) = coding_fr;
                inp(inp <= fr_thre/1000) = 1;
            end

            inp(inp < 1) = 0; 

            for i = 1:inp_neuron_num
                if inp(i) == 1 
                    inp_spike_flag(i) = 1;
                    inp_spike_count(i) = inp_spike_count(i) + 1;
                    inp_spike_time(i,inp_spike_count(i)) = t;
                end
                if inp_spike_flag(i) == 1
                    if inp_spike_count(i) < 2
                        inp_res_ker(i,rem(t,inp_res_ker_length) + 1) =  exp(-(t-inp_spike_time(i,inp_spike_count(i)))/tau_f);
                    else
                        inp_last_spike(i) = t - inp_spike_time(i,inp_spike_count(i) - 1);
                        inp_res_ker(i,rem(t,inp_res_ker_length) + 1) =  exp(-inp_last_spike(i)/tau_f) + exp(-(t-inp_spike_time(i,inp_spike_count(i)))/tau_f);
                    end
                end
            end

            % E neurons
            exc_u(:,1) = (inp_res_ker(:,rem(t-5,inp_res_ker_length) + 1)' * inp_exc_weight + exc_res_ker(:,rem(t-5,exc_res_ker_length) + 1)' * exc_exc_updated_w(:,:,rem(t-5,exc_exc_length) + 1)...
                 - inh_res_ker(:,rem(t-2,inh_res_ker_length) + 1)' * inh_exc_weight)'; 

            exc_ref_temp = max(min(round(gamrnd(ref_shape,exc_ref_mu,exc_neuron_num,1)),exc_ref_max),1); 
            exc(:,1:end-1) = exc(:,2:end); 
            for j = 1:exc_neuron_num
                if exc_u(j,1) > u_thre && max(exc(j,end-exc_ref_temp(j)+1:end)) == 0 
                    exc(j,end) = 1; 
                else
                    exc(j,end) = 0;
                end

                if exc(j,end) == 1 
                    exc_spike_flag(j) = 1;
                    exc_spike_count(j) = exc_spike_count(j) + 1;
                    exc_spike_time(j,exc_spike_count(j)) = t;
                end

                if exc_spike_flag(j) == 1
                    if exc_spike_count(j) < 2
                        exc_res_ker(j,rem(t,exc_res_ker_length) + 1) = exp(-(t-exc_spike_time(j,exc_spike_count(j)))/tau_f);
                    else
                        exc_last_spike(j) = t - exc_spike_time(j,exc_spike_count(j) - 1);
                        exc_res_ker(j,rem(t,exc_res_ker_length) + 1) = exp(-exc_last_spike(j)/tau_f) + exp(-(t-exc_spike_time(j,exc_spike_count(j)))/tau_f);
                    end
                end
            end    

            if strcmp(plasticity,'ON')
                % E->E STDP
                r1_exc_exc(:,:,rem(t,2) + 1) = exc_exc_weight_mask .* (r1_exc_exc(:,:,rem(t-1,2) + 1) + (r1_exc_exc_0 - r1_exc_exc(:,:,rem(t-1,2) + 1))/tau_r1);
                o1_exc_exc(:,:,rem(t,2) + 1) = exc_exc_weight_mask .* (o1_exc_exc(:,:,rem(t-1,2) + 1) + (o1_exc_exc_0 - o1_exc_exc(:,:,rem(t-1,2) + 1))/tau_o1);
                r2_exc_exc(:,:,rem(t,2) + 1) = exc_exc_weight_mask .* (r2_exc_exc(:,:,rem(t-1,2) + 1) + (r2_exc_exc_0 - r2_exc_exc(:,:,rem(t-1,2) + 1))/tau_r2);
                o2_exc_exc(:,:,rem(t,2) + 1) = exc_exc_weight_mask .* (o2_exc_exc(:,:,rem(t-1,2) + 1) + (o2_exc_exc_0 - o2_exc_exc(:,:,rem(t-1,2) + 1))/tau_o2);   
                for j = 1:exc_neuron_num
                    for i = 1:exc_neuron_num
                        if exc_exc_weight(i,j) > 0
                            if exc_spike_count(i) ~= 0 && t == exc_spike_time(i,exc_spike_count(i))
                                r1_exc_exc(i,j,rem(t,2) + 1) = r1_exc_exc(i,j,rem(t-1,2) + 1) + 1;
                                r2_exc_exc(i,j,rem(t,2) + 1) = r2_exc_exc(i,j,rem(t-1,2) + 1) + 1;
                                exc_exc_updated_w(i,j,rem(t,exc_exc_length) + 1) = -o1_exc_exc(i,j,rem(t,2) + 1)*(A2_minus_exc + A3_minus_exc*r2_exc_exc(i,j,rem(t-1,2) + 1)) + exc_exc_updated_w(i,j,rem(t-1,exc_exc_length) + 1);
                            elseif exc_spike_count(j) ~= 0 && t == exc_spike_time(j,exc_spike_count(j))
                                o1_exc_exc(i,j,rem(t,2) + 1) = o1_exc_exc(i,j,rem(t-1,2) + 1) + 1;
                                o2_exc_exc(i,j,rem(t,2) + 1) = o2_exc_exc(i,j,rem(t-1,2) + 1) + 1;
                                exc_exc_updated_w(i,j,rem(t,exc_exc_length) + 1) = r1_exc_exc(i,j,rem(t,2) + 1)*(A2_plus_exc + A3_plus_exc*o2_exc_exc(i,j,rem(t-1,2) + 1)) + exc_exc_updated_w(i,j,rem(t-1,exc_exc_length) + 1);
                            else
                                exc_exc_updated_w(i,j,rem(t,exc_exc_length) + 1) = exc_exc_updated_w(i,j,rem(t-1,exc_exc_length) + 1);
                            end

                            if exc_exc_updated_w(i,j,rem(t,exc_exc_length) + 1) > max_exc_exc_weight
                                exc_exc_updated_w(i,j,rem(t,exc_exc_length) + 1) = max_exc_exc_weight; 
                            elseif exc_exc_updated_w(i,j,rem(t,exc_exc_length) + 1) < min_exc_exc_weight 
                                exc_exc_updated_w(i,j,rem(t,exc_exc_length) + 1) = min_exc_exc_weight;
                            end
                        end
                    end
                end

                if rem(t,20) == 0 % Homeostatic plasticity
                    exc_exc_updated_w_temp = exc_exc_updated_w(:,:,rem(t,exc_exc_length) + 1);
                    exc_exc_updated_w_temp = exc_exc_updated_w_temp - (sum(exc_exc_updated_w_temp) - exc_exc_weight_0)./exc_exc_weight_pre_num;
                    exc_exc_updated_w_temp(exc_exc_updated_w_temp > max_exc_exc_weight) = max_exc_exc_weight;
                    exc_exc_updated_w_temp(exc_exc_updated_w_temp < min_exc_exc_weight) = min_exc_exc_weight;
                    exc_exc_updated_w_temp = exc_exc_updated_w_temp .* exc_exc_weight_mask;
                    exc_exc_updated_w(:,:,rem(t,exc_exc_length) + 1) = exc_exc_updated_w_temp;
                end
            end

            % I neurons
            inh_u(:,1) = (inp_res_ker(:,rem(t-2,inp_res_ker_length) + 1)' * inp_inh_weight + exc_res_ker(:,rem(t-2,exc_res_ker_length) + 1)' * exc_inh_weight - inh_res_ker(:,rem(t-2,inh_res_ker_length) + 1)' * inh_inh_weight)'; % 加入synaptic delay

            inh_ref_temp = max(min(round(gamrnd(ref_shape,inh_ref_mu,inh_neuron_num,1)),inh_ref_max),1);
            inh(:,1:end-1) = inh(:,2:end);
            for k = 1:inh_neuron_num         
                if inh_u(k,1) > u_thre && max(inh(k,end-inh_ref_temp(k)+1:end)) == 0
                    inh(k,end) = 1; 
                else
                    inh(k,end) = 0;
                end

                if inh(k,end) == 1 
                    inh_spike_flag(k) = 1;
                    inh_spike_count(k) = inh_spike_count(k) + 1;
                    inh_spike_time(k,inh_spike_count(k)) = t;
                end

                if inh_spike_flag(k) == 1
                    if inh_spike_count(k) < 2
                        inh_res_ker(k,rem(t,inh_res_ker_length) + 1) =  exp(-(t-inh_spike_time(k,inh_spike_count(k)))/tau_f);
                    else
                        inh_last_spike(k) = t - inh_spike_time(k,inh_spike_count(k) - 1);
                        inh_res_ker(k,rem(t,inh_res_ker_length) + 1) =  exp(-inh_last_spike(k)/tau_f) + exp(-(t-inh_spike_time(k,inh_spike_count(k)))/tau_f);
                    end
                end
            end


            if rem(t,Nskip) == 0
                exc_exc_weight_save(:,:,floor(t/Nskip) + 1) = exc_exc_updated_w(:,:,rem(t,exc_exc_length) + 1); 
            end

        end
    end


    % Plot and save firing rate pattern
    figure,
    imagesc(inp_fr),colorbar()
    xlabel('t (s)')
    ylabel('Neuron index')
    title('Input neuron');
    saveas(gcf,[Save_figure_dir,'AI_inp=',num2str(inp_neuron_num),'_',num2str(single_pattern_num),'P_',num2str(single_pattern_show_time),'S_',...
        num2str(asso_pattern_show_time),'A_inp_fr_inp-E=',num2str(inp_exc_ini),...
        '_E-E=',num2str(exc_exc_ini),'_E-I=',num2str(exc_inh_ini),'_I-E=',num2str(inh_exc_ini),'_I-I=',num2str(inh_inh_ini),'_pattern_time=',num2str(pattern_time),'_interval_time=',num2str(interval_time),...
        '_LTP_lr=',num2str(A2_plus_exc_novel),'_LTD_lr=',num2str(A2_minus_exc_novel),'_exc_ref=',num2str(exc_ref_mu),'_trail=',num2str(trail),'.jpg']);  %保存当前窗口的图像

    figure,
    imagesc(exc_fr),colorbar()
    xlabel('t (s)')
    ylabel('Neuron index')
    title('Exc neuron');
    saveas(gcf,[Save_figure_dir,'AI_inp=',num2str(inp_neuron_num),'_',num2str(single_pattern_num),'P_',num2str(single_pattern_show_time),'S_',...
        num2str(asso_pattern_show_time),'A_exc_fr_inp-E=',num2str(inp_exc_ini),...
        '_E-E=',num2str(exc_exc_ini),'_E-I=',num2str(exc_inh_ini),'_I-E=',num2str(inh_exc_ini),'_I-I=',num2str(inh_inh_ini),'_pattern_time=',num2str(pattern_time),'_interval_time=',num2str(interval_time),...
        '_LTP_lr=',num2str(A2_plus_exc_novel),'_LTD_lr=',num2str(A2_minus_exc_novel),'_exc_ref=',num2str(exc_ref_mu),'_trail=',num2str(trail),'.jpg']);  %保存当前窗口的图像

    figure,
    imagesc(inh_fr),colorbar()
    xlabel('t (s)')
    ylabel('Neuron index')
    title('Inh neuron');
    saveas(gcf,[Save_figure_dir,'AI_inp=',num2str(inp_neuron_num),'_',num2str(single_pattern_num),'P_',num2str(single_pattern_show_time),'S_',...
        num2str(asso_pattern_show_time),'A_inh_fr_inp-E=',num2str(inp_exc_ini),...
        '_E-E=',num2str(exc_exc_ini),'_E-I=',num2str(exc_inh_ini),'_I-E=',num2str(inh_exc_ini),'_I-I=',num2str(inh_inh_ini),'_pattern_time=',num2str(pattern_time),'_interval_time=',num2str(interval_time),...
        '_LTP_lr=',num2str(A2_plus_exc_novel),'_LTD_lr=',num2str(A2_minus_exc_novel),'_exc_ref=',num2str(exc_ref_mu),'_trail=',num2str(trail),'.jpg']);  %保存当前窗口的图像

    % Save firing rate and weight file
    save([Save_data_dir,'AI_inp=',num2str(inp_neuron_num),'_',num2str(single_pattern_num),'P_',num2str(single_pattern_show_time),'S_',...
        num2str(asso_pattern_show_time),'A_inp_fr_inp-E=',num2str(inp_exc_ini),...
        '_E-E=',num2str(exc_exc_ini),'_E-I=',num2str(exc_inh_ini),'_I-E=',num2str(inh_exc_ini),'_I-I=',num2str(inh_inh_ini),'_pattern_time=',num2str(pattern_time),'_interval_time=',num2str(interval_time),...
        '_LTP_lr=',num2str(A2_plus_exc_novel),'_LTD_lr=',num2str(A2_minus_exc_novel),'_exc_ref=',num2str(exc_ref_mu),'_trail=',num2str(trail),'.mat'],'inp_fr');
    save([Save_data_dir,'AI_inp=',num2str(inp_neuron_num),'_',num2str(single_pattern_num),'P_',num2str(single_pattern_show_time),'S_',...
        num2str(asso_pattern_show_time),'A_exc_fr_inp-E=',num2str(inp_exc_ini),...
        '_E-E=',num2str(exc_exc_ini),'_E-I=',num2str(exc_inh_ini),'_I-E=',num2str(inh_exc_ini),'_I-I=',num2str(inh_inh_ini),'_pattern_time=',num2str(pattern_time),'_interval_time=',num2str(interval_time),...
        '_LTP_lr=',num2str(A2_plus_exc_novel),'_LTD_lr=',num2str(A2_minus_exc_novel),'_exc_ref=',num2str(exc_ref_mu),'_trail=',num2str(trail),'.mat'],'exc_fr');
    save([Save_data_dir,'AI_inp=',num2str(inp_neuron_num),'_',num2str(single_pattern_num),'P_',num2str(single_pattern_show_time),'S_',...
        num2str(asso_pattern_show_time),'A_inh_fr_inp-E=',num2str(inp_exc_ini),...
        '_E-E=',num2str(exc_exc_ini),'_E-I=',num2str(exc_inh_ini),'_I-E=',num2str(inh_exc_ini),'_I-I=',num2str(inh_inh_ini),'_pattern_time=',num2str(pattern_time),'_interval_time=',num2str(interval_time),...
        '_LTP_lr=',num2str(A2_plus_exc_novel),'_LTD_lr=',num2str(A2_minus_exc_novel),'_exc_ref=',num2str(exc_ref_mu),'_trail=',num2str(trail),'.mat'],'inh_fr');
    save([Save_data_dir,'AI_inp=',num2str(inp_neuron_num),'_',num2str(single_pattern_num),'P_',num2str(single_pattern_show_time),'S_',...
        num2str(asso_pattern_show_time),'A_EE_weight_inp-E=',num2str(inp_exc_ini),...
        '_E-E=',num2str(exc_exc_ini),'_E-I=',num2str(exc_inh_ini),'_I-E=',num2str(inh_exc_ini),'_I-I=',num2str(inh_inh_ini),'_pattern_time=',num2str(pattern_time),'_interval_time=',num2str(interval_time),...
        '_LTP_lr=',num2str(A2_plus_exc_novel),'_LTD_lr=',num2str(A2_minus_exc_novel),'_exc_ref=',num2str(exc_ref_mu),'_trail=',num2str(trail),'.mat'],'exc_exc_weight_save');

    % Save the spike time of excitatory neurons
    save([Save_data_dir,'AI_inp=',num2str(inp_neuron_num),'_',num2str(single_pattern_num),'P_',num2str(single_pattern_show_time),'S_',...
        num2str(asso_pattern_show_time),'A_exc_spike_time_inp-E=',num2str(inp_exc_ini),...
        '_E-E=',num2str(exc_exc_ini),'_E-I=',num2str(exc_inh_ini),'_I-E=',num2str(inh_exc_ini),'_I-I=',num2str(inh_inh_ini),'_pattern_time=',num2str(pattern_time),'_interval_time=',num2str(interval_time),...
        '_LTP_lr=',num2str(A2_plus_exc_novel),'_LTD_lr=',num2str(A2_minus_exc_novel),'_exc_ref=',num2str(exc_ref_mu),'_trail=',num2str(trail),'.mat'],'exc_spike_time');

    print_info = '总运行时间%s\n';
    fprintf(print_info,(datetime('now') - start_time)) % print run time

    close all
end





