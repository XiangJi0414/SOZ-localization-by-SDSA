%% read the original data
clear;
subname = 'S09';
report_label = {'A 01';'A 02';'A 03';'D 01';'D 02';'D 03';'E 01';'E 02';'E 03'};% label was provided by clinicians
filename1 = 'S09_data_epi.mat';
load(filename1);
cfg = [];
cfg.viewmode = 'vertical';
cfg.channel = 'all';
ft_databrowser(cfg, data_epi); % view the SEEG data

%% compute SDSA matrix
sample = data_epi.fsample;
label = data_epi.label;
win_length = 1; % 1s一个时间窗口，计算一秒内的标注差
Num_segment = size(data_epi.trial,2);
if Num_segment == 1
    start_time = data_epi.time{1,1}(1);
    end_time = data_epi.time{1,1}(end);
    time_label = start_time:1:(end_time-1);
    std_jibo = []; % SDSA matrix
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

end

%% locate the SOZ for S09 based on SDSA
if Num_segment == 1

    ph_jibo = smooth(mean_jibo,30, 'sgolay');
    gd_jibo = gradient(ph_jibo);
    [~,onset_idx] = max(gd_jibo);
    onset_time = time_label(onset_idx);
    disp('The T_b for seizure boosting:');
    disp(onset_time);

    csum_jibo = sum(std_jibo(:,1:onset_idx),2);
    [bingzao_ds,bingzao_idx] = sort(csum_jibo,'descend');

    disp(bingzao_label);

    class2 = kmeans(bingzao_ds,3);
    class1 = bingzao_label(find(class2 == class2(1)));
    disp('SOZ by kmeans：');
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
        fprintf('The T_b for No.%d seizure boosting: %d',j,onset_time(1,j));

        csum_jibo = sum(std_jibo{1,j}(:,1: onset_idx),2);
        [bingzao_ds{1,j},bingzao_idx{1,j}] = sort(csum_jibo,'descend');
%         figure(Num_segment+4+j);
%         bar(bingzao_ds{1,j});
        bingzao_label{1,j} = label(bingzao_idx{1,j});
%         xticks(1:length(label));
%         xticklabels(bingzao_label{1,j});
        disp(bingzao_label{1,j});
        class2 = kmeans(bingzao_ds{1,j},3);
        class1{1,j} = bingzao_label{1,j}(find(class2 == class2(1)));
        disp('SOZ by kmeans：')
        disp(class1{1,j});
        if j == 1
            classduolei = class1{1,j};
        else
            classduolei = union(classduolei,class1{1,j});
        end
    end
end

%% plot figures
probility = [];
for j = 1:1:Num_segment
    max_j = max(bingzao_ds{1,j});
    probility(:,j) = bingzao_ds{1,j}/max_j;
    figure(2*Num_segment+4+j);
    bar(probility(:,j));
    ax = gca;
    xticks(1:length(bingzao_label{1,j}));
    xticklabels(bingzao_label{1,j});

    ax.XAxis.FontSize = 8;
    %     xticklabels.FontSize = 3;
    ylabel('SOZ probility','FontSize',18);
    xlabel('Channels','FontSize',18);
    title('SOZ probability of each channel','FontSize',18);
end

%% figure 1D
j = 1;
figure(14);
plot(time_label{1,j},mean_jibo{1,j},'Color','#A30234','LineStyle','-');
hold on;
plot(time_label{1,j},ph_jibo{1,j},'Color','#00545f','LineStyle','-','LineWidth',1);
%     plot(onset_time(1,j),ph_jibo{1,j}(time_label{1,j}==onset_time(1,j)),'Color','#7E2F8E','Marker','v','MarkerFaceColor','#7E2F8E','MarkerSize',8);
plot(onset_time(1,j),ph_jibo{1,j}(time_label{1,j}==onset_time(1,j)),'Color','r','Marker','v','MarkerFaceColor','r','MarkerSize',8);
ax = gca;
ylim([0,250]);
xticks([15538,15614,15690,15766,15842]);
xticklabels([15538,15614,15690,15766,15842]);
ax.XAxis.FontSize = 10;
yticks([50,100,150,200,250]);
yticklabels([50,100,150,200,250]);
ax.YAxis.FontSize = 10;
title('Average SDSA over all channels during the Seizure 1','FontSize',18,'FontName','Calibri');
xlabel('time/s','FontName','Calibri','FontSize',15);
ylabel('SDSA_a_v_g','FontName','Calibri','FontSize',15);
legend('raw SDSA_a_v_g','smoothing SDSA_a_v_g','T_b','FontSize',12,'FontName','Calibri');
set(gca,'ygrid','on');
%% figure 1E
j =1;
figure(15);
imagesc(std_jibo{1,j});
colorbar;
ax = gca;
set(gca,'YTickLabel',[]);
xlabel('time/s','FontName','Calibri','FontSize',15);
set(gca,'XTick',[],'XTickLabel',[]);
hold on;
% xline(find(time_label{1,j}==onset_time(1,j)),'Color','r','LineStyle','--','LineWidth',1);
ylabel('Channels','FontName','Calibri','FontSize',15);
title('SDSA of each channel during the Seizure 1','FontName','Calibri','FontSize',18);
%% Figure 1F
figure(16);
j = 1;
[~,onset_idx] = max(gd_jibo{1,j});
csum_jibo = sum(std_jibo{1,j}(:,1: onset_idx),2);
barh(csum_jibo,'EdgeColor','k','FaceColor','#0076C0');
ax=gca;
%         bingzao_label{1,j} = label(bingzao_idx{1,j});
set(gca,'Ydir','reverse');
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
ylabel('Channels','FontName','Calibri','FontSize',15);
xlabel('SDSA Integration','FontName','Calibri','FontSize',15);
title('SOZ probability of each channel','FontName','Calibri','FontSize',18);
%% Figure 1G
AUC_idx = ismember(bingzao_label{1,1},report_label);
perf_label = bingzao_label{1,1};
perf_label(AUC_idx) = {'tar'};
[perf_X,perf_Y,perf_T,perf_AUC] = perfcurve(perf_label,bingzao_ds{1,1},'tar');
disp('AUC is：');
disp(perf_AUC);
figure(17);
plot(perf_X,perf_Y,'Color','#E37C1D','LineWidth',1.5);
xlabel('False Positive Rate','FontName','Calibri','FontSize',15);
ylabel('True Positive Rate','FontName','Calibri','FontSize',15);
title('Accuracy of SOZ probability rank, AUC=0.9980','FontName','Calibri','FontSize',18);hold on;
area(perf_X,perf_Y,'FaceColor',[227 124 29]/255,'FaceAlpha',0.1);hold off;
%% Figure 1H
figure(18);
h = bar(1:10,probility(1:10,1));

hold on;
h1 = bar(11:19,probility(11:19,1));
h2 = bar(20:118,probility(20:118,1));

ax = gca;
set(gca,'XTick',[]);
set(gca,'XTickLabel',[]);
%     xticklabels.FontSize = 3;
ylabel('SOZ probability','FontName','Calibri','FontSize',15);
xlabel('Channels sorted by descend SOZ probability','FontName','Calibri','FontSize',15);

for j = 1 : 10
    h(1, 1).FaceColor = 'flat';
    h(1, 1).CData(j,:) = [206,128,128]/255;
    h(1, 1).EdgeColor = 'k';
end
for j = 1:9
    h1(1, 1).FaceColor = 'flat';
    h1(1, 1).CData(j,:) = [171,180,125]/255;
    h1(1, 1).EdgeColor = 'k';
end
for j = 1:99
    h2(1, 1).FaceColor = 'flat';
    h2(1, 1).CData(j,:) = [161,197,203]/255;
    h2(1,1).EdgeColor = 'k';
end
title('Classify channels by k-means clustering','FontName','Calibri','FontSize',18);