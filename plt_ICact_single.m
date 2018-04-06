subjectId = 1;
% The components here are using index 
IC_ori_selected = [8,9,11,12,13,14,15,17,18,19,22];
IC_clean_selected = [2,4,5,8,3,1,10,14,18,16,15,11];
thresASR = [1, 2.5, 3, 5, 10, 20, 30, 40, 50, 70, 100, 200, 500, 1000];

% load parameters
addpath('parameters/');
load(sprintf('paras_s%02d.mat',subjectId));
eval(sprintf('paras = paras_s%02d;',subjectId));
nIC = size(paras.ori_Ch_var.ori_var_wholeSec,1);

saveFigPath = sprintf('figure/s%02d/',subjectId);
if ~exist(saveFigPath,'dir')
    mkdir(saveFigPath);
end

act_ICclean_wholeSec = zeros(nIC,length(thresASR));
act_ICori_wholeSec = zeros(nIC,length(thresASR));
act_ICclean_cleanSec = zeros(nIC,length(thresASR));
act_ICori_cleanSec = zeros(nIC,length(thresASR));

for i = 1:length(thresASR)
    act_ICclean_wholeSec(:,i) = paras.ASR_result(i).whole_Section.ICact_ori_clean;
    act_ICori_wholeSec(:,i) = paras.ASR_result(i).whole_Section.ICact_ori;
    act_ICclean_cleanSec(:,i) = paras.ASR_result(i).clean_Section.ICact_ori_clean;
    act_ICori_cleanSec(:,i) = paras.ASR_result(i).clean_Section.ICact_ori;
end

%% plot
lw = 3;
linestyle = {'-','--',':','-.','-'};
ori_ICori_wholeSec = paras.ori_ICA_result.whole_Section.ICact_ori;
ori_ICclean_wholeSec = paras.ori_ICA_result.whole_Section.ICact_ori_clean;
ori_ICori_cleanSec = paras.ori_ICA_result.clean_Section.ICact_ori;
ori_ICclean_cleanSec = paras.ori_ICA_result.clean_Section.ICact_ori_clean;

%% plot all IC_ori in whole section (raw)
cmap = jet(nIC);
figure()
hold on
for i = 1:nIC
    plot(thresASR, act_ICori_wholeSec(i,:), 'color', cmap(i,:))
end
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('ASR cutoff parameter')
ylabel('Power of source activities(\muV^2)');
title('IC_{ori} in wholeSec (raw)');
grid on

%% plot all IC_ori_clean in whole section (raw)
cmap = jet(nIC);
figure()
hold on
for i = 1:nIC
    plot(thresASR, act_ICclean_wholeSec(i,:), 'color', cmap(i,:))
end
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('ASR cutoff parameter')
ylabel('Power of source activities(\muV^2)');
title('IC_{ori,clean} in wholeSec (raw)');
grid on

%% plot all IC_ori in clean section (raw)
cmap = jet(nIC);
figure()
hold on
for i = 1:nIC
    plot(thresASR, act_ICori_cleanSec(i,:), 'color', cmap(i,:))
end
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('ASR cutoff parameter')
ylabel('Power of source activities(\muV^2)');
title('IC_{ori} in wholeSec (raw)');
grid on

%% plot all IC_ori_clean in clean section (raw)
cmap = jet(nIC);
figure()
hold on
for i = 1:nIC
    plot(thresASR, act_ICclean_cleanSec(i,:), 'color', cmap(i,:))
end
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('ASR cutoff parameter')
ylabel('Power of source activities(\muV^2)');
title('IC_{ori,clean} in wholeSec (raw)');
grid on

%% plot all IC_ori in whole section (%)
cmap = jet(nIC);
act = act_ICori_wholeSec./ori_ICori_wholeSec*100;
figure()
hold on
for i = 1:nIC
    plot(thresASR, act(i,:), 'color', cmap(i,:))
end
set(gca,'xscale','log')
xlabel('ASR cutoff parameter')
ylabel('Retained power(%)');
title('IC_{ori} in wholeSec (%)');
grid on

%% plot all IC_ori_clean in whole section (%)
cmap = jet(nIC);
act = act_ICclean_wholeSec./ori_ICclean_wholeSec*100;
figure()
hold on
for i = 1:nIC
    plot(thresASR, act(i,:), 'color', cmap(i,:))
end
set(gca,'xscale','log')
xlabel('ASR cutoff parameter');
ylabel('Retained power(%)');
title('IC_{ori,clean} in wholeSec (%)');
grid on

%% plot all IC_ori in clean section (%)
cmap = jet(nIC);
act = act_ICori_cleanSec./ori_ICori_cleanSec*100;
% act = act_ICori_cleanSec./act_ICori_cleanSec(:,14)*100;
figure()
hold on
set(gca,'xscale','log')
xlabel('ASR cutoff parameter');
ylabel('Retained power(%)');
title('IC_{ori} in cleanSec (%)');
grid on
for i = 1:nIC
    plot(thresASR, act(i,:), 'color', cmap(i,:))
end


%% plot all IC_ori_clean in clean section (%)
cmap = jet(nIC);
act = act_ICclean_cleanSec./ori_ICclean_cleanSec*100;
figure()
hold on
for i = 1:nIC
    plot(thresASR, act(i,:), 'color', cmap(i,:))
end
set(gca,'xscale','log')
xlabel('ASR cutoff parameter');
ylabel('Retained power(%)');
title('IC_{ori,clean} in cleanSec (%)');
grid on



%% IC_ori
cmap = jet(length(IC_ori_selected));
figure()
hold on
for i = 1:length(IC_ori_selected)
    Hb = plot(thresASR,act_ICori_wholeSec(IC_ori_selected(i),:),'DisplayName',sprintf('comp%d',IC_ori_selected(i)),'color',cmap(i,:));
%     set(Hb,'LineStyle',linestyle{i});
    set(Hb,'LineWidth',lw);
end

set(gca,'xscale','log')
set(gca,'fontsize',20);    
xlabel('ASR cutoff parameter')
ylabel('Power of source activities(\muV^2)');
title('IC_{ori} in wholeSec (raw)');
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gca,'yscale','log');
legend(findobj(gca,'-regexp','DisplayName', '[^'']'));
grid on

%% IC_ori_clean
cmap = jet(length(IC_clean_selected));
figure()
hold on
for i = 1:length(IC_clean_selected)
    Hb = plot(thresASR,act_ICclean_wholeSec(IC_clean_selected(i),:),'DisplayName',sprintf('comp%d',IC_clean_selected(i)),'color',cmap(i,:));
%     set(Hb,'LineStyle',linestyle{i});
    set(Hb,'LineWidth',lw);
end

set(gca,'xscale','log')
set(gca,'fontsize',20);    
xlabel('ASR cutoff parameter')
ylabel('Power of source activities(\muV^2)');
title('IC_{ori,clean} in wholeSec (raw)');
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gca,'yscale','log');
legend(findobj(gca,'-regexp','DisplayName', '[^'']'));
grid on


