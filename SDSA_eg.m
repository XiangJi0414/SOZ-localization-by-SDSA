% % % This script tells how to compute SDSA matrix from SEEG data and then compute SDSA_avg vector
% % % smoothen the SDSA_avg vector to find T_b
% % % integration from T_s to T_b to locate SOZ

%% Read the original ictal SEEG data, here We use random functions to generate raw data. 
% Assuming 118 channels, F_sample 512Hz, start_time 0s, end_time 600s. 

clear;
data_epi.trial{1,1} = rand([118,600*512]); % Replace to your real SEEG data in practice.
sample = 512;
label = cell(118,1);
win_length = 1; % one window per second.
Num_segment = size(data_epi.trial,2);
%% Compute SDSA matrix
if Num_segment == 1
    start_time = 0;
    end_time = 600;
    time_label = start_time:1:(end_time-1);
    std_jibo = [];
    for i = 1:1:end_time-start_time
        window = data_epi.trial{1, 1}(:,(i-1)*win_length*sample+1:i*win_length*sample);
        std_dev = std(window,0,2);
        std_jibo(:,i) = std_dev;
    end
    mean_jibo = mean(std_jibo,1);
    figure(2);
    plot(time_label,mean_jibo,'r-');
    title('Average SDSA over all channels during the seizure');
    xlabel('time/s','FontSize',13);
    ylabel('Standard deviation of spike amplitude (SDSA)','FontSize',13);
    figure(3);
    imagesc(std_jibo);
    colorbar;
    ax = gca;
    yticks(1:length(label));
    yticklabels(label);
    ax.TickLabelInterpreter = 'none';
    ax.YAxis.FontSize = 4;
    xlabel('time/s','FontSize',13);
    ylabel('Channels','FontSize',13);
    title('Standard deviation of spike amplitude (SDSA) of each channel during the seizure','FontSize',13);
else
    std_jibo = cell(1,Num_segment);
    mean_jibo = cell(1,Num_segment);
    start_time = zeros(1,Num_segment);
    end_time = zeros(1,Num_segment);
    time_label = cell(1,Num_segment);
    Time_freq = cell(1,Num_segment);
    for j = 1:1:Num_segment
        start_time(1,j) = data_epi.time{1,j}(1);
        end_time(1,j) = data_epi.time{1,j}(end);
        time_label{1,j} = start_time(1,j):1:(end_time(1,j)-1);
        std_jibo{1,j} = [];
        for i = 1:1:end_time(1,j)-start_time(1,j)
            window = data_epi.trial{1, j}(:,(i-1)*win_length*sample+1:i*win_length*sample);
            std_dev = std(window,0,2);
            std_jibo{1,j}(:,i) = std_dev;
        end
        mean_jibo{1,j} = mean(std_jibo{1,j},1);
    end
    %%% plot
    figure(2);
    for j = 1:1:Num_segment
        subplot(2,2,j);
        plot(time_label{1,j},mean_jibo{1,j},'r-');
        title('Average SDSA over all channels during the seizure','FontSize',13);
        xlabel('time/s','FontSize',13);
        ylabel('Standard deviation of spike amplitude (SDSA)','FontSize',13);

    end
    figure(3);
    for j = 1:1:Num_segment
        subplot(1,2,j);
        imagesc(std_jibo{1,j});
        colorbar;
        ax = gca;
        yticks(1:length(label));
        yticklabels(label);
        ax.TickLabelInterpreter = 'none';
        ax.YAxis.FontSize = 4;
        xlabel('time','FontSize',13);

        ylabel('Channels','FontSize',13);
        title('Standard deviation of spike amplitude (SDSA) of each channel during the seizure','FontSize',13);
    end
end
%% smoothen the SDSA_avg vector to find T_b
if Num_segment == 1

    ph_jibo = smooth(mean_jibo,30, 'sgolay');
    gd_jibo = gradient(ph_jibo);
    [~,onset_idx] = max(gd_jibo);
    onset_time = time_label(onset_idx);
    disp('T_b:');
    disp(onset_time);
%   Integration from Ts to Tb
    csum_jibo = sum(std_jibo(:,1:onset_idx),2);
    [bingzao_ds,bingzao_idx] = sort(csum_jibo,'descend');
    figure(6);
    bar(bingzao_ds);
    bingzao_label = label(bingzao_idx);
    xticks(1:length(label));
    xticklabels(bingzao_label);
    disp(bingzao_label); % descend sorting
    % clustering 
    class2 = kmeans(bingzao_ds,3);
    class1 = bingzao_label(find(class2 == class2(1)));
    disp('kmeans cluster for SOZ:');
    disp(class1);
else
    ph_jibo = cell(1,Num_segment);
    gd_jibo = cell(1,Num_segment);
    onset_time = zeros(1,Num_segment);
    bingzao_ds = cell(1,Num_segment);
    bingzao_idx = cell(1,Num_segment);
    bingzao_label = cell(1,Num_segment);
    class1 =  cell(1,Num_segment);
    for j = 1:1:Num_segment

        ph_jibo{1,j} = smooth(mean_jibo{1,j},30,'sgolay');
        gd_jibo{1,j} = gradient(ph_jibo{1,j});
        [~,onset_idx] = max(gd_jibo{1,j});
        onset_time(1,j) = time_label{1,j}(onset_idx);
        fprintf('No. %d T_b:%d',j,onset_time(1,j));

        csum_jibo = sum(std_jibo{1,j}(:,1: onset_idx),2);
        [bingzao_ds{1,j},bingzao_idx{1,j}] = sort(csum_jibo,'descend');
        figure(Num_segment+4+j);
        bar(bingzao_ds{1,j});
        bingzao_label{1,j} = label(bingzao_idx{1,j});
        xticks(1:length(label));
        xticklabels(bingzao_label{1,j});
        disp(bingzao_label{1,j});
        class2 = kmeans(bingzao_ds{1,j},3);
        class1{1,j} = bingzao_label{1,j}(find(class2 == class2(1)));
        disp('kmeans cluster for SOZ:')
        disp(class1{1,j});
        if j == 1
            classduolei = class1{1,j};
        else
            classduolei = union(classduolei,class1{1,j});
        end
    end
end
