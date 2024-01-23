close all

clearvars -except Behav_path Path_figures Path_variables Code_path MEG_path

% plot Figure 3 - 2

SubjectID={'S01' 'S02' 'S04' 'S05' 'S06' 'S07' 'S08' 'S10' 'S11' 'S12' 'S15' 'S17' 'S18' 'S19' 'S20' 'S22' 'S24','S21','S25','S23'};
Drug=[3 1 1 2 3 2 3; 3 1 3 2 1 2 0; 1 3 2 1 3 2 0; 2 1 1 3 2 3 0; 2 1 3 1 2 3 0; 1 3 2 3 1 2 0; 3 1 1 2 2 3 0; 1 2 2 3 1 3 1; 1 2 2 3 3 1 1; 3 1 2 1 3 2 0; 1 3 3 2 2 1 0; 2 1 3 2 3 1 2; 1 2 3 2 1 3 1; 2 3 1 2 1 3 0; 1 2 3 1 3 2 0; 3 1 2 2 1 3 0; 3 1 1 2 3 2 0; 2 1 1 2 3 3 0; 3 2 1 3 1 2 0;2 3 2 1 1 3 0];
Stimulus=[40 36.75 36.75 40 40 40 40; 40 36.75 36.75 40 40 40 40; 40 36.75 36.75 36.75 40 40 40; 40 36.75 36.75 36.75 40 40 40; 36.75 40 40 40 40 40 40; 40 36.75 36.75 36.75 40 40 40; 40 36.75 36.75 36.75 40 40 40; 40 36.75 36.75 40 40 40 40; 40 36.75 36.75 36.75 40 40 40; 40 36.75 36.75 36.75 40 40 40; 40 36.75 36.75 40 40 40 40; 40 36.75 36.75 36.75 40 40 40; 40 36.75 36.75 36.75 40 40 40; 40 40 40 40 40 40 40; 40 40 40 40 40 40 40; 40 40 40 40 40 40 40; ...
    40 40 40 40 40 40 40;40 40 40 40 40 40 40;40 40 40 40 40 40 40;40 40 40 40 40 40 40];
Freq=[29:0.25:48];

addpath(Code_path)

%% Total Power

load([Path_variables,'/Figure3/GrandAvg_drug_A_EAC_ITPC.mat'])
load([Path_variables,'/Figure3/GrandAvg_drug_B_EAC_ITPC.mat'])
load([Path_variables,'/Figure3/GrandAvg_drug_C_EAC_ITPC.mat'])

%% cluster permutation test for TFR of ITPC between drugs

ITC_drug_A.powspctrm=permute(ITC_drug_A.powspctrm,[4 1 2 3]);
ITC_drug_B.powspctrm=permute(ITC_drug_B.powspctrm,[4 1 2 3]);
ITC_drug_C.powspctrm=permute(ITC_drug_C.powspctrm,[4 1 2 3]);

TFR_drug_A=ITC_drug_A;
TFR_drug_B=ITC_drug_B;
TFR_drug_C=ITC_drug_C;


TFR_A = (mean(ITC_drug_A.powspctrm,1));
TFR_B = (mean(ITC_drug_B.powspctrm,1));
TFR_C = (mean(ITC_drug_C.powspctrm,1));

TFR_drug_AB     = TFR_A-TFR_B;
TFR_drug_AB     = reshape(TFR_drug_AB,size(TFR_drug_A.freq,2),size(TFR_drug_A.time,2));

TFR_drug_AC     = TFR_A-TFR_C;
TFR_drug_AC     = reshape(TFR_drug_AC,size(TFR_drug_A.freq,2),size(TFR_drug_A.time,2));

TFR_drug_CB     = TFR_C-TFR_B;
TFR_drug_CB     = reshape(TFR_drug_CB,size(TFR_drug_A.freq,2),size(TFR_drug_A.time,2));

cfg = [];
cfg.channel          = 'all';
cfg.latency          = [0.2 2.225];  
cfg.frequency        = [-5 5];
cfg.avgoverchan      = 'yes';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.tail             = 0;
cfg.correctm         = 'cluster';
cfg.clusterstatistic = 'maxsum';
cfg.clustertail      = 0;
cfg.clusteralpha     = 0.05;
cfg.correcttail      = 'prob';
cfg.alpha            = 0.05;
cfg.numrandomization = 10000;

subj = size(SubjectID,2);
design = zeros(2,2*subj);
for i = 1:subj
  design(1,i) = i;
end
for i = 1:subj
  design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;

[statAB] = ft_freqstatistics(cfg, TFR_drug_A, TFR_drug_B);
save([Path_variables,'/Figure3/EAC_TFR_DrugAB_ITPC.mat'],'statAB')
[statAC] = ft_freqstatistics(cfg, TFR_drug_A, TFR_drug_C);
save([Path_variables,'/Figure3/EAC_TFR_DrugAC_ITPC.mat'],'statAC')
[statCB] = ft_freqstatistics(cfg, TFR_drug_C, TFR_drug_B);
save([Path_variables,'/Figure3/EAC_TFR_DrugCB_ITPC.mat'],'statCB')

load([Path_variables,'/Figure3/EAC_TFR_DrugAB_ITPC.mat'])
load([Path_variables,'/Figure3/EAC_TFR_DrugAC_ITPC.mat'])
load([Path_variables,'/Figure3/EAC_TFR_DrugCB_ITPC.mat'])


pvalueAB=reshape(statAB.mask,size(statAB.freq,2),size(statAB.time,2));
pvalueAB=double(pvalueAB);
maskAB=pvalueAB;

pvalueAC=reshape(statAC.mask,size(statAC.freq,2),size(statAC.time,2));
pvalueAC=double(pvalueAC);
maskAC=pvalueAC;

pvalueCB=reshape(statCB.mask,size(statCB.freq,2),size(statCB.time,2));
pvalueCB=double(pvalueCB);
maskCB=pvalueCB;

figure
%sgtitle('Total Power')
s1 = subplot(3,1,1);
imagesc(TFR_drug_A.time,TFR_drug_A.freq,TFR_drug_AB) 
xline(0,'--k')
xline(2,'--k')
hold on
[~,c] = contour(statAB.time,statAB.freq,maskAB,[1 1],'k');
c.LineWidth = 1;
xlim([-0.25 2.225])
set(gca,'xtick',[])
set(gca,'xticklabel',[])
ylim([-5 5])
yticks([-5 0 5])
caxis([-0.15 0.15])
c = colorbar;
title('Lorazepam - Placebo')
set(gca,'YDir','normal')
box off
pbaspect([1.7 1 1])
ax=gca;
ax.FontSize = 7;

s2 = subplot(3,1,2);
imagesc(TFR_drug_A.time,TFR_drug_A.freq,TFR_drug_AC) 
xline(0,'--k')
xline(2,'--k')
hold on
[~,c] = contour(statAB.time,statAB.freq,maskAC,[1 1],'k');
c.LineWidth = 1;
xlim([-0.25 2.225])
set(gca,'xtick',[])
set(gca,'xticklabel',[])
ylim([-5 5])
yticks([-5 0 5])
ylabel('Hz from StimFrequency')
caxis([-0.15 0.15])
c = colorbar;
c.Label.String = '% Power change';
c.Label.FontSize = 7;
title('Lorazepam - Memantine')
set(gca,'YDir','normal')
box off
pbaspect([1.7 1 1])
ax=gca;
ax.FontSize = 7;

s3 = subplot(3,1,3);
imagesc(TFR_drug_A.time,TFR_drug_A.freq,TFR_drug_CB) 
xline(0,'--k')
xline(2,'--k')
hold on
[M,c] = contour(statAB.time,statAB.freq,maskCB,[1 1],'k');
c.LineWidth = 1;
xlabel('Time from stimulus onset (s)')
xlim([-0.25 2.225])
xticks([0 1 2])
ylim([-5 5])
caxis([-0.15 0.15])
c = colorbar;
title('Memantine - Placebo')
set(gca,'YDir','normal')
box off
pbaspect([1.7 1 1])
ax=gca;
ax.FontSize = 7;
%set(gcf,'units','points','position',[910,421,540,350])

set(s1,'Units','centimeters', 'position', [2.56822222222222 9.25787459150327 4.4686 2.6288]);
set(s2,'Units','centimeters', 'position', [2.56822222222222 5.346840073529413 4.4686 2.6288]);
set(s3,'Units','centimeters', 'position', [2.56822222222222 1.435805555555556 4.4686 2.6288]);

cd(Path_figures)

saveas(gcf,'Figure3_2_B.png')


%% RoT

freq=[-7:0.25:7];
time=[-0.2:0.025:(2.5-0.0025)];

cfg                  = [];
cfg.channel          = 'all';
cfg.latency          = [0.2 2.225];
cfg.frequency        = [-1 1];
cfg.avgoverchan      = 'yes';
cfg.avgoverfreq      = 'yes';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';  
cfg.tail             = 0; 
cfg.correctm         = 'cluster';
cfg.clusterstatistic = 'maxsum';  
cfg.clustertail      = 0;
cfg.clusteralpha     = 0.05;  
cfg.correcttail      = 'prob';
cfg.alpha            = 0.05;
cfg.numrandomization = 10000;

subj = size(SubjectID,2);
design = zeros(2,2*subj);
for i = 1:subj
  design(1,i) = i;
end
for i = 1:subj
  design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;

[stat_AB] = ft_freqstatistics(cfg, TFR_drug_A, TFR_drug_B);
save([Path_variables,'/Figure3/EAC_TFR_DrugAB_ITPC_RoT.mat'],'statAB')
[stat_CA] = ft_freqstatistics(cfg, TFR_drug_A, TFR_drug_C);
save([Path_variables,'/Figure3/EAC_TFR_DrugAC_ITPC_RoT.mat'],'statCA')
[stat_CB] = ft_freqstatistics(cfg, TFR_drug_C, TFR_drug_B);
save([Path_variables,'/Figure3/EAC_TFR_DrugCB_ITPC_RoT.mat'],'statCB')

load([Path_variables,'/Figure3/EAC_TFR_DrugAB_ITPC_RoT.mat'])
load([Path_variables,'/Figure3/EAC_TFR_DrugAC_ITPC_RoT.mat'])
load([Path_variables,'/Figure3/EAC_TFR_DrugCB_ITPC_RoT.mat'])




id_freqs = [find(freq==cfg.frequency(1)):find(freq==cfg.frequency(2))];
id_time  = [29:110];

stat_AB.d_mask=double(stat_AB.mask);
stat_AB.d_mask=reshape(stat_AB.d_mask,1,82);
stat_CB.d_mask=double(stat_CB.mask);
stat_CB.d_mask=reshape(stat_CB.d_mask,1,82);
stat_CA.d_mask=double(stat_CA.mask);
stat_CA.d_mask=reshape(stat_CA.d_mask,1,82);

stat_AB.d_mask(stat_AB.d_mask==0)=nan;
stat_CB.d_mask(stat_CB.d_mask==0)=nan;
stat_CA.d_mask(stat_CA.d_mask==0)=nan;

stat_AB.mask_ok=nan(1,120);
stat_AB.mask_ok(1,id_time)=stat_AB.d_mask;
stat_CB.mask_ok=nan(1,120);
stat_CB.mask_ok(1,id_time)=stat_CB.d_mask;
stat_CA.mask_ok=nan(1,120);
stat_CA.mask_ok(1,id_time)=stat_CA.d_mask;

TFR__A=(nanmean(nanmean(TFR_drug_A.powspctrm(:,:,25:33,:),2),3));
TFRA=reshape(TFR__A,size(SubjectID,2),size(TFR_drug_A.time,2));

TFR__B=(nanmean(nanmean(TFR_drug_B.powspctrm(:,:,25:33,:),2),3));
TFRB=reshape(TFR__B,size(SubjectID,2),size(TFR_drug_A.time,2));

TFR__C=(nanmean(nanmean(TFR_drug_C.powspctrm(:,:,25:33,:),2),3));
TFRC=reshape(TFR__C,size(SubjectID,2),size(TFR_drug_A.time,2));

TFRhann_A=nanmean(nanmean(nanmean(TFR_drug_A.powspctrm(:,:,25:33,:),2),3),1);
TFRhann_A=reshape(TFRhann_A,1,size(TFR_drug_C.time,2));

TFRhann_B=nanmean(nanmean(nanmean(TFR_drug_B.powspctrm(:,:,25:33,:),2),3),1);
TFRhann_B=reshape(TFRhann_B,1,size(TFR_drug_C.time,2));

TFRhann_C=nanmean(nanmean(nanmean(TFR_drug_C.powspctrm(:,:,25:33,:),2),3),1);
TFRhann_C=reshape(TFRhann_C,1,size(TFR_drug_C.time,2));

sz=30;

s1=figure;
scatter(1:120,stat_AB.mask_ok*40+8,sz,[0 162 255]./255,'filled')
hold on
scatter(1:120,stat_CA.mask_ok*40+6,sz,[255/2 255/2 255/2]./255,'filled')
hold on
scatter(1:120,stat_CB.mask_ok*40+4,sz,[238 34 12]./255,'filled')
hold on
plot(TFR_drug_A.time, TFRhann_A,'Color',[0 162 255]./255)
shadedErrorBar(1:size(TFRhann_A,2),TFRhann_A,nanstd(TFRA)./sqrt(size(SubjectID,2)),'lineProps',{'Color',[0 162 255]./255});
hold on
plot(TFR_drug_A.time, TFRhann_B,'Color',[255/2 255/2 255/2]./255)
shadedErrorBar(1:size(TFRhann_B,2),TFRhann_B,nanstd(TFRB)./sqrt(size(SubjectID,2)),'lineProps',{'Color',[255/2 255/2 255/2]./255});
hold on
plot(TFR_drug_A.time, TFRhann_C,'Color',[238 34 12]./255)
shadedErrorBar(1:size(TFRhann_C,2),TFRhann_C,nanstd(TFRC)./sqrt(size(SubjectID,2)),'lineProps',{'Color',[238 34 12]./255});
hold on
xline(21,'--k')
xline(101,'--k')
yline(0,'--k')
xlabel('Time from stimulus onset (s)')
xticks([21 61 101])
xticklabels({'0','1','2'})
xlim([11 110])
ylabel('% Power change')
ylim([0 0.40]);
yticks([0 0.15 0.30])
legend('','','','Lorazepam (LO)','','Placebo (PL)','','Memantine (ME)','','','','Location','south')
pbaspect([1.5 1 1])
box off
ax=gca;
ax.FontSize = 7;
%set(gcf,'units','points','position',[455,421,340,200])
set(s1,'Units','centimeters', 'position', [1.55927777777778 0.776111111111111 8.3873 5.5915]);

cd(Path_figures)

saveas(gcf,'Figure3_2_C.png')

%% Violin transient

cfg                  = [];
cfg.channel          = 'all';
cfg.latency          = [0.1 0.4];
cfg.frequency        = [-1 1];
cfg.avgoverchan      = 'yes';
cfg.avgovertime      = 'yes';
cfg.avgoverfreq      = 'yes';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT'; 
cfg.tail             = 0; 
cfg.correctm         = 'no';
%cfg.clusterstatistic = 'maxsum';
%cfg.clustertail      = 0;
%cfg.clusteralpha     = 0.05;
cfg.correcttail      = 'prob';
cfg.alpha            = 0.05;
cfg.numrandomization = 10000;

subj = size(SubjectID,2);
design = zeros(2,2*subj);
for i = 1:subj
  design(1,i) = i;
end
for i = 1:subj
  design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;

[stat_AB] = ft_freqstatistics(cfg, TFR_drug_A, TFR_drug_B);
save([Path_variables,'/Figure3/EAC_TFR_DrugAB_ITPC_Violin1.mat'],'statAB')

pvalueAB  = num2str(stat_AB.prob);
[stat_CB] = ft_freqstatistics(cfg, TFR_drug_C, TFR_drug_B);
save([Path_variables,'/Figure3/EAC_TFR_DrugCB_ITPC_Violin1.mat'],'statCB')

pvalueCB  = num2str(stat_CB.prob);
[stat_CA] = ft_freqstatistics(cfg, TFR_drug_A, TFR_drug_C);
save([Path_variables,'/Figure3/EAC_TFR_DrugCA_ITPC_Violin1.mat'],'statCA')

load([Path_variables,'/Figure3/EAC_TFR_DrugAB_ITPC_Violin1.mat'])
load([Path_variables,'/Figure3/EAC_TFR_DrugCB_ITPC_Violin1.mat'])
load([Path_variables,'/Figure3/EAC_TFR_DrugCA_ITPC_Violin1.mat'])



pvalueAC  = num2str(stat_CA.prob);

TFR_A_pow = permute(TFR_drug_A.powspctrm,[3,4,1,2]);
TFR_A_pow = TFR_A_pow;
TFR_B_pow = permute(TFR_drug_B.powspctrm,[3,4,1,2]);
TFR_B_pow = TFR_B_pow;
TFR_C_pow = permute(TFR_drug_C.powspctrm,[3,4,1,2]);
TFR_C_pow = TFR_C_pow;

Subjects(1,:) = reshape(mean(mean(TFR_B_pow(25:33,25:37,:),2),1),size(SubjectID,2),1);
Subjects(2,:) = reshape(mean(mean(TFR_A_pow(25:33,25:37,:),2),1),size(SubjectID,2),1);
Subjects(3,:) = reshape(mean(mean(TFR_C_pow(25:33,25:37,:),2),1),size(SubjectID,2),1);

colors = [[255/2 255/2 255/2]./255;[0 162 255]./255;[238 34 12]./255];

figure
violin([reshape(mean(mean(TFR_B_pow(25:33,25:37,:),2),1),size(SubjectID,2),1),...
   reshape(mean(mean(TFR_A_pow(25:33,25:37,:),2),1),size(SubjectID,2),1),...
   reshape(mean(mean(TFR_C_pow(25:33,25:37,:),2),1),size(SubjectID,2),1)],'facecolor',colors,'edgecolor','k','mc','k','medc','k--','L')
hold on
for sub=1:size(SubjectID,2)
    for dr=1:3
        scatter(dr,Subjects(dr,sub),7,'k','filled')
        hold on
    end
end
text(1.5,100,pvalueAB)
text(2,120,pvalueCB)
text(2.5,100,pvalueAC)
%ylabel('% Power change','FontSize', 11)
yline(0,'--k')
%xticks([1,2,3])
ylim([0 0.60])
yticks([0 0.30 0.60])
%xticklabels({'Lorazepam','Placebo','Memantine'})
set(gca,'xtick',[])
set(gca,'xticklabel',[])
title('Transient Component')
box off
ax=gca;
ax.FontSize = 7;
set(gcf,'units','points','position',[455,421,150,140])


cd(Path_figures)

saveas(gcf,'Figure3_2_D1.png')


%% Violin sustained

cfg                  = [];
cfg.channel          = 'all';
cfg.latency          = [0.65 1.75];
cfg.frequency        = [-1 1];
cfg.avgoverchan      = 'yes';
cfg.avgovertime      = 'yes';
cfg.avgoverfreq      = 'yes';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT'; 
cfg.tail             = 0; 
cfg.correctm         = 'no';
%cfg.clusterstatistic = 'maxsum';
%cfg.clustertail      = 0;
%cfg.clusteralpha     = 0.05;
cfg.correcttail      = 'prob';
cfg.alpha            = 0.05;
cfg.numrandomization = 10000;

subj = size(SubjectID,2);
design = zeros(2,2*subj);
for i = 1:subj
  design(1,i) = i;
end
for i = 1:subj
  design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;

[stat_AB] = ft_freqstatistics(cfg, TFR_drug_A, TFR_drug_B);
save([Path_variables,'/Figure3/EAC_TFR_DrugAB_ITPC_Violin2.mat'],'statAB')
pvalueAB  = num2str(stat_AB.prob);
[stat_CB] = ft_freqstatistics(cfg, TFR_drug_C, TFR_drug_B);
save([Path_variables,'/Figure3/EAC_TFR_DrugCB_ITPC_Violin2.mat'],'statCB')
pvalueCB  = num2str(stat_CB.prob);
[stat_CA] = ft_freqstatistics(cfg, TFR_drug_A, TFR_drug_C);
save([Path_variables,'/Figure3/EAC_TFR_DrugCA_ITPC_Violin2.mat'],'statCA')
pvalueAC  = num2str(stat_CA.prob);


load([Path_variables,'/Figure3/EAC_TFR_DrugAB_ITPC_Violin2.mat'])
load([Path_variables,'/Figure3/EAC_TFR_DrugCB_ITPC_Violin2.mat'])
load([Path_variables,'/Figure3/EAC_TFR_DrugCA_ITPC_Violin2.mat'])


TFR_A_pow = permute(TFR_drug_A.powspctrm,[3,4,1,2]);
TFR_A_pow = TFR_A_pow;
TFR_B_pow = permute(TFR_drug_B.powspctrm,[3,4,1,2]);
TFR_B_pow = TFR_B_pow;
TFR_C_pow = permute(TFR_drug_C.powspctrm,[3,4,1,2]);
TFR_C_pow = TFR_C_pow;

Subjects(1,:) = reshape(mean(mean(TFR_B_pow(25:33,47:91,:),2),1),size(SubjectID,2),1);
Subjects(2,:) = reshape(mean(mean(TFR_A_pow(25:33,47:91,:),2),1),size(SubjectID,2),1);
Subjects(3,:) = reshape(mean(mean(TFR_C_pow(25:33,47:91,:),2),1),size(SubjectID,2),1);

colors = [[255/2 255/2 255/2]./255;[0 162 255]./255;[238 34 12]./255];

figure
violin([reshape(mean(mean(TFR_B_pow(25:33,47:91,:),2),1),size(SubjectID,2),1),...
   reshape(mean(mean(TFR_A_pow(25:33,47:91,:),2),1),size(SubjectID,2),1),...
   reshape(mean(mean(TFR_C_pow(25:33,47:91,:),2),1),size(SubjectID,2),1)],'facecolor',colors,'edgecolor','k','mc','k','medc','k--','L')
hold on
for sub=1:size(SubjectID,2)
    for dr=1:3
        scatter(dr,Subjects(dr,sub),7,'k','filled')
        hold on
    end
end
text(1.5,100,pvalueAB)
text(2,120,pvalueCB)
text(2.5,100,pvalueAC)
yline(0,'--k')
%xticks([1,2,3])
ylim([0 0.60])
yticks([0 0.30 0.60])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
set(gca,'xtick',[])
set(gca,'xticklabel',[])
%ylabel('% Power change','FontSize', 11)
legend('','Mean','Median','Location','southeast')
%xticklabels({'Lorazepam','Placebo','Memantine'})
title('Sustained Component')
box off
ax=gca;
ax.FontSize = 7;
set(gcf,'units','points','position',[455,421,150,140])
cd(Path_figures)

saveas(gcf,'Figure3_2_D2.png')

