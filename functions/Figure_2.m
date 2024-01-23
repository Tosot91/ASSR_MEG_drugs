%% plot figure 2

clearvars -except Behav_path Path_figures Path_variables

close all

SubjectID={'S01' 'S02' 'S04' 'S05' 'S06' 'S07' 'S08' 'S10' 'S11' 'S12' 'S15' 'S17' 'S18' 'S19' 'S20' 'S22' 'S24','S21','S25','S23'};
Drug=[3 1 1 2 3 2 3; 3 1 3 2 1 2 0; 1 3 2 1 3 2 0; 2 1 1 3 2 3 0; 2 1 3 1 2 3 0; 1 3 2 3 1 2 0; 3 1 1 2 2 3 0; 1 2 2 3 1 3 1; 1 2 2 3 3 1 1; 3 1 2 1 3 2 0; 1 3 3 2 2 1 0; 2 1 3 2 3 1 2; 1 2 3 2 1 3 1; 2 3 1 2 1 3 0; 1 2 3 1 3 2 0; 3 1 2 2 1 3 0; 3 1 1 2 3 2 0; 2 1 1 2 3 3 0; 3 2 1 3 1 2 0;2 3 2 1 1 3 0];
Stimulus=[40 36.75 36.75 40 40 40 40; 40 36.75 36.75 40 40 40 40; 40 36.75 36.75 36.75 40 40 40; 40 36.75 36.75 36.75 40 40 40; 36.75 40 40 40 40 40 40; 40 36.75 36.75 36.75 40 40 40; 40 36.75 36.75 36.75 40 40 40; 40 36.75 36.75 40 40 40 40; 40 36.75 36.75 36.75 40 40 40; 40 36.75 36.75 36.75 40 40 40; 40 36.75 36.75 40 40 40 40; 40 36.75 36.75 36.75 40 40 40; 40 36.75 36.75 36.75 40 40 40; 40 40 40 40 40 40 40; 40 40 40 40 40 40 40; 40 40 40 40 40 40 40; ...
    40 40 40 40 40 40 40;40 40 40 40 40 40 40;40 40 40 40 40 40 40;40 40 40 40 40 40 40];
Freq=[29:0.25:48];

load([Path_variables,'/Figure2/Significant_Channels_stat.mat'])

% get significant sensors

Significant_channels=[];
if isfield(stat,'posclusterslabelmat')
    for i=1:size(stat.posclusters,2)
        if stat.posclusters(i).prob<0.05
            Significant_channels=unique([Significant_channels;stat.label(find(stat.posclusterslabelmat==i))]);
        end
    end
end


for i=1:size(Significant_channels,1)
    
    string=cell2mat(Significant_channels(i));
    string((end-2):end)=[];
    Significant_channels_ok{i,:}=string;
    
end
Significant_channels=unique(Significant_channels_ok);

%% plot TFR at sensor level for placebo sessions

load([Path_variables,'/Figure2/40Hz_PL_Sensor_all.mat'])
TFR_40_PL_all           = TFR_drug_B;
TFR_40_PL_all.powspctrm = (nanmean(TFR_40_PL_all.powspctrm,4))*100;
load([Path_variables,'/Figure2/40Hz_PL_Sensor_sig.mat'])
TFR_40_PL_sig           = TFR_drug_B;
TFR_40_PL_sig.powspctrm = (nanmean(TFR_40_PL_sig.powspctrm,4))*100;
TFR_40                  = TFR_40_PL_sig.powspctrm;
TFR_40_sig              = mean(TFR_40,1);
TFR_40_sig              = reshape(TFR_40_sig,size(TFR_40_PL_sig.freq,2),size(TFR_40_PL_sig.time,2));

clims = [-30 30];

figure
s1 = subplot(1,5,1);
imagesc(TFR_40_PL_sig.time,TFR_40_PL_sig.freq,TFR_40_sig,clims)
xline(0,'--k')
xline(2,'--k')
yline(40,'--k')
yline(36.75,'--k')
xlabel('Time (s)')
xlim([-0.25 2.225])
xticks([0 1 2])
caxis([-30 30])
c = colorbar('location','eastoutside');
c.Label.String = '% Power change';
c.Label.FontSize = 7;
ylabel('Frequency (Hz)')
yticks([5 36.75 40 60])
title('Total Power')
set(gca,'YDir','normal')
box off
pbaspect([1 2 1])
ax=gca;
ax.FontSize = 7;

set(s1,'Units','centimeters', 'position', [2.71514378962627 1.62983333333333 2.0113 4.0217]);
%set(gcf,'units','points','position',[910,421,460,320])

cd(Path_figures)
saveas(gcf,'Figure2_A_2.png')

f=figure;
s1 = subplot(3,1,1);
cfg                     = [];
cfg.channel             = 'all';
cfg.markersymbol        = '.';
cfg.markercolor         = [1 1 1];
cfg.highlight           = 'on';
cfg.highlightchannel    = Significant_channels;
cfg.highlightsymbol     = '.';
cfg.highlightcolor      = [0 0 0];
cfg.highlightsize       = 7;
cfg.zlim                = [-30 30];
cfg.xlim                = [0.1 2];
cfg.ylim                = [39 41];
cfg.layout              = 'CTF275_helmet';
cfg.colormap            = 'jet';
%cfg.colorbar            = 'EastOutside';
cfg.marker              = 'on';
cfg.title               = 'Placebo';
cfg.comment             = 'no';
cfg.interactive         = 'no';
ft_topoplotTFR(cfg, TFR_40_PL_all);
f.Position(3:4) = [500 360];

cd(Path_figures)
saveas(gcf,'Figure2_A.png')

%% EAC Placebo

load([Path_variables,'/Figure2/40Hz_PL_EAC.mat'])
TFRhann_drug40_EAC=TFRhann_drug;
TFRhann_drug40_EAC.powspctrm=(mean(TFRhann_drug40_EAC.powspctrm,4))*100;
TFRhann_drug40_EAC.powspctrm=reshape(TFRhann_drug40_EAC.powspctrm,size(TFRhann_drug.freq,2),size(TFRhann_drug.time,2));
load([Path_variables,'/Figure2/40Hz_PL_phase_EAC.mat'])
TFRhann_drug40_EAC_phase=TFRhann_drug_phase;
TFRhann_drug40_EAC_phase.powspctrm=(mean(TFRhann_drug40_EAC_phase.powspctrm,4))*100;
TFRhann_drug40_EAC_phase.powspctrm=reshape(TFRhann_drug40_EAC_phase.powspctrm,size(TFRhann_drug.freq,2),size(TFRhann_drug.time,2));
load([Path_variables,'/Figure2/37Hz_PL_EAC.mat'])
TFRhann_drug37_EAC=TFRhann_drug;
TFRhann_drug37_EAC.powspctrm=(mean(TFRhann_drug37_EAC.powspctrm,4))*100;
TFRhann_drug37_EAC.powspctrm=reshape(TFRhann_drug37_EAC.powspctrm,size(TFRhann_drug.freq,2),size(TFRhann_drug.time,2));
load([Path_variables,'/Figure2/37Hz_PL_phase_EAC.mat'])
TFRhann_drug37_EAC_phase=TFRhann_drug_phase;
TFRhann_drug37_EAC_phase.powspctrm=(mean(TFRhann_drug37_EAC_phase.powspctrm,4))*100;
TFRhann_drug37_EAC_phase.powspctrm=reshape(TFRhann_drug37_EAC_phase.powspctrm,size(TFRhann_drug.freq,2),size(TFRhann_drug.time,2));
load([Path_variables,'/Figure2/40Hz_PL_dlPFC.mat'])
TFRhann_drug40_dlPFC=TFRhann_drug;
TFRhann_drug40_dlPFC.powspctrm=(mean(TFRhann_drug40_dlPFC.powspctrm,4))*100;
TFRhann_drug40_dlPFC.powspctrm=reshape(TFRhann_drug40_dlPFC.powspctrm,size(TFRhann_drug.freq,2),size(TFRhann_drug.time,2));

clims = [-30 30];

figure
s1 = subplot(1,5,1);
imagesc(TFRhann_drug40_EAC.time,TFRhann_drug40_EAC.freq,TFRhann_drug40_EAC.powspctrm,clims)
xline(0,'--k')
xline(2,'--k')
yline(36.75,'--k')
yline(40,'--k')
xlim([-0.25 2.225])
xticks([0 1 2])
caxis([-30 30])
ylabel('Frequency (Hz)')
yticks([5 36.75 40 60])
title('Total Power')
set(gca,'YDir','normal')
box off
pbaspect([1 2 1])
ax=gca;
ax.FontSize = 7;

s2 = subplot(1,5,2);
imagesc(TFRhann_drug40_EAC.time,TFRhann_drug40_EAC.freq,TFRhann_drug37_EAC.powspctrm,clims)
xline(0,'--k')
xline(2,'--k')
yline(36.75,'--k')
yline(40,'--k')
xlim([-0.25 2.225])
set(gca,'XDir','normal')
xticks([0 1 2])
caxis([-30 30])
title('Total Power')
set(gca,'YDir','normal')
yticks([5 36.75 40 60])
box off
pbaspect([1 2 1])
ax=gca;
ax.FontSize = 7;

s3 = subplot(1,5,3);
imagesc(TFRhann_drug40_EAC.time,TFRhann_drug40_EAC.freq,TFRhann_drug40_dlPFC.powspctrm,clims)
xline(0,'--k')
xline(2,'--k')
yline(36.75,'--k')
yline(40,'--k')
xlim([-0.25 2.225])
xticks([0 1 2])
caxis([-30 30])
yticks([5 36.75 40 60])
xlabel('Time from stimulus onset (s)','FontSize', 11)
title('Total Power','FontSize', 11)
set(gca,'YDir','normal')
box off
pbaspect([1 2 1])
ax=gca;
ax.FontSize = 7;

s4 = subplot(1,5,4);
imagesc(TFRhann_drug40_EAC.time,TFRhann_drug40_EAC.freq,TFRhann_drug40_EAC_phase.powspctrm,clims)
xline(0,'--k')
xline(2,'--k')
yline(36.75,'--k')
yline(40,'--k')
xlim([-0.25 2.225])
set(gca,'XDir','normal')
xticks([0 1 2])
caxis([-30 30])
title('Phase-Locked')
set(gca,'YDir','normal')
yticks([5 36.75 40 60])
box off
pbaspect([1 2 1])
ax=gca;
ax.FontSize = 7;

s5 = subplot(1,5,5);
imagesc(TFRhann_drug40_EAC.time,TFRhann_drug40_EAC.freq,TFRhann_drug37_EAC_phase.powspctrm,clims)
xline(0,'--k')
xline(2,'--k')
yline(36.75,'--k')
yline(40,'--k')
xlim([-0.25 2.225])
xticks([0 1 2])
caxis([-30 30])
c = colorbar('location','eastoutside');
c.Label.String = '% Power change';
c.Label.FontSize = 7;
yticks([5 36.75 40 60])
title('Phase-Locked')
set(gca,'YDir','normal')
box off
pbaspect([1 2 1])
ax=gca;
ax.FontSize = 7;

set(s1,'Units','centimeters', 'position', [2.71514378962627 1.62983333333333 2.0113 4.0217]);
set(s2,'Units','centimeters', 'position', [5.93164705760013 1.62983333333333 2.0113 4.0217]);
set(s3,'Units','centimeters', 'position', [9.148150325573987 1.62983333333333 2.0113 4.0217]);
set(s4,'Units','centimeters', 'position', [12.364653593547844 1.62983333333333 2.0113 4.0217]);
set(s5,'Units','centimeters', 'position', [15.434235294117649 1.62983333333333 2.0113 4.0217]);

cd(Path_figures)
saveas(gcf,'Figure2_C_E.png')