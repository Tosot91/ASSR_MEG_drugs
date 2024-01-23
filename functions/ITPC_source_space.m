close all

clearvars -except Behav_path Path_figures Path_variables Code_path MEG_path

SubjectID={'S01' 'S02' 'S05' 'S15' 'S19' 'S20' 'S22' 'S24','S21','S25','S23'};
Drug=[3 1 1 2 3 2 3; 3 1 3 2 1 2 0; 2 1 1 3 2 3 0; 1 3 3 2 2 1 0; 2 3 1 2 1 3 0; 1 2 3 1 3 2 0; 3 1 2 2 1 3 0; 3 1 1 2 3 2 0; 2 1 1 2 3 3 0; 3 2 1 3 1 2 0;2 3 2 1 1 3 0];
Stimulus=[40 36.75 36.75 40 40 40 40; 40 36.75 36.75 40 40 40 40; 40 36.75 36.75 36.75 40 40 40; 40 36.75 36.75 40 40 40 40; 40 40 40 40 40 40 40; 40 40 40 40 40 40 40; 40 40 40 40 40 40 40; 40 40 40 40 40 40 40;40 40 40 40 40 40 40;40 40 40 40 40 40 40;40 40 40 40 40 40 40];
Freq=[5:0.25:60];

%% calculate ITPC in source space for PL sessions w/ 40Hz AM stimulus

for Type=2
  
    TFR_this_drug=[];
    TFR_this_drug_phase=[];
     
    for sub=1:size(SubjectID,2)
        Sessions = find(Drug(sub,:)==Type);
        StimSessions = Stimulus(sub,Sessions);

        TFR_this_subject=[];
        TFR_this_subject_phase=[];

        if isempty(Sessions)==0
            
            for Session=1:size(Sessions,2)
                if Session==3
                 load([MEG_path,'/source/',cell2mat(SubjectID(sub)),'/SES11/A1_ERF_l_long.mat'])
                    datas1=comb_dict{1,1};
                    load([MEG_path,'/source/',cell2mat(SubjectID(sub)),'/SES11/LBelt_ERF_l_long.mat'])
                    datas2=comb_dict{1,1};
                    load([MEG_path,'/source/',cell2mat(SubjectID(sub)),'/SES11/MBelt_ERF_l_long.mat'])
                    datas3=comb_dict{1,1};
                    load([MEG_path,'/source/',cell2mat(SubjectID(sub)),'/SES11/PBelt_ERF_l_long.mat'])
                    datas4=comb_dict{1,1};
                    load([MEG_path,'/source/',cell2mat(SubjectID(sub)),'/SES11/RI_ERF_l_long.mat'])
                    datas5=comb_dict{1,1};
                    load([MEG_path,'/source/',cell2mat(SubjectID(sub)),'/SES11/A1_ERF_r_long.mat'])
                    datas6=comb_dict{1,1};
                    load([MEG_path,'/source/',cell2mat(SubjectID(sub)),'/SES11/LBelt_ERF_r_long.mat'])
                    datas7=comb_dict{1,1};
                    load([MEG_path,'/source/',cell2mat(SubjectID(sub)),'/SES11/MBelt_ERF_r_long.mat'])
                    datas8=comb_dict{1,1};
                    load([MEG_path,'/source/',cell2mat(SubjectID(sub)),'/SES11/PBelt_ERF_r_long.mat'])
                    datas9=comb_dict{1,1};
                    load([MEG_path,'/source/',cell2mat(SubjectID(sub)),'/SES11/RI_ERF_r_long.mat'])
                    datas10=comb_dict{1,1};
                else
                   load([MEG_path,'/source/',cell2mat(SubjectID(sub)),'/SES',num2str(Sessions(Session)),'/A1_ERF_l_long.mat'])
                    datas1=comb_dict{1,1};
                    load([MEG_path,'/source/',cell2mat(SubjectID(sub)),'/SES',num2str(Sessions(Session)),'/LBelt_ERF_l_long.mat'])
                    datas2=comb_dict{1,1};
                    load([MEG_path,'/source/',cell2mat(SubjectID(sub)),'/SES',num2str(Sessions(Session)),'/MBelt_ERF_l_long.mat'])
                    datas3=comb_dict{1,1};
                    load([MEG_path,'/source/',cell2mat(SubjectID(sub)),'/SES',num2str(Sessions(Session)),'/PBelt_ERF_l_long.mat'])
                    datas4=comb_dict{1,1};
                    load([MEG_path,'/source/',cell2mat(SubjectID(sub)),'/SES',num2str(Sessions(Session)),'/RI_ERF_l_long.mat'])
                    datas5=comb_dict{1,1};
                    load([MEG_path,'/source/',cell2mat(SubjectID(sub)),'/SES',num2str(Sessions(Session)),'/A1_ERF_r_long.mat'])
                    datas6=comb_dict{1,1};
                    load([MEG_path,'/source/',cell2mat(SubjectID(sub)),'/SES',num2str(Sessions(Session)),'/LBelt_ERF_r_long.mat'])
                    datas7=comb_dict{1,1};
                    load([MEG_path,'/source/',cell2mat(SubjectID(sub)),'/SES',num2str(Sessions(Session)),'/MBelt_ERF_r_long.mat'])
                    datas8=comb_dict{1,1};
                    load([MEG_path,'/source/',cell2mat(SubjectID(sub)),'/SES',num2str(Sessions(Session)),'/PBelt_ERF_r_long.mat'])
                    datas9=comb_dict{1,1};
                    load([MEG_path,'/source/',cell2mat(SubjectID(sub)),'/SES',num2str(Sessions(Session)),'/RI_ERF_r_long.mat'])
                    datas10=comb_dict{1,1};
                end

                datas_all = cat(2,datas1.erfdata,datas2.erfdata,datas3.erfdata,datas4.erfdata,datas5.erfdata,datas6.erfdata,datas7.erfdata,datas8.erfdata,datas9.erfdata,datas10.erfdata);

               if Session==3
                    load([MEG_path,'/sensors/',cell2mat(SubjectID(sub)),'/SES11/MNE_data_clean_postICA_',cell2mat(SubjectID(sub)),'SES11.mat'])
                else
                    load([MEG_path,'/sensors/',cell2mat(SubjectID(sub)),'/SES',num2str(Sessions(Session)),'/MNE_data_clean_postICA_',cell2mat(SubjectID(sub)),'SES',num2str(Sessions(Session)),'.mat'])
                end
                
                data.trial = datas_all;
                data.label = data.label(1:size(datas_all,2),:);
    
                StimFreq = StimSessions(1,Session);

                % take sessions with stimulus 40 Hz
                if StimFreq==40

                    id_ripple = find(data.trialinfo(:,1)==1);

                    cfg              = [];
                    cfg.toilim       = [-0.5 2.5-0.0025];
                    data             = ft_redefinetrial(cfg,data);

                    cfg              = [];
                    cfg.trials       = id_ripple;
                    ERF              = ft_timelockanalysis(cfg,data);

                    cfg              = [];
                    cfg.output       = 'fourier';
                    cfg.channel      = 'MEG';
                    cfg.trials       = id_ripple;
                    cfg.keeptrials   = 'yes';
                    cfg.method       = 'mtmconvol';
                    cfg.taper        = 'hanning';
                    cfg.foi          = Freq;     
                    cfg.pad          = 4;
                    cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   
                    cfg.toi          = -0.5:0.025:(2.5-0.0025);            
                    freq = ft_freqanalysis(cfg, data);

                    cfg              = [];
                    cfg.output       = 'pow';
                    cfg.channel      = 'MEG';
                    cfg.trials       = id_ripple;
                    cfg.keeptrials   = 'yes';
                    cfg.method       = 'mtmconvol';
                    cfg.taper        = 'hanning';
                    cfg.foi          = Freq;     
                    cfg.pad          = 4;
                    cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   
                    cfg.toi          = -0.5:0.025:(2.5-0.0025);            
                    TFR_pow = ft_freqanalysis(cfg, data);

                    
                
                itc           = [];
                itc.label     = freq.label;
                itc.freq      = freq.freq;
                itc.time      = freq.time;
                itc.dimord    = 'chan_freq_time';
    
                F = freq.fourierspctrm;   
                N = size(F,1);          
    
                % compute inter-trial phase coherence (itpc)
                itc.itpc      = F./abs(F);         % divide by amplitude
                itc.itpc      = sum(itc.itpc,1);   % sum angles
                itc.itpc      = abs(itc.itpc)/N;   % take the absolute value and normalize
                itc.itpc      = squeeze(itc.itpc); % remove the first singleton dimension

                itc.itpc=mean(itc.itpc,1);
                itc.label={'EAC'};

                    
                    % concatenate ITPCs
                    Relch = itc.itpc;
                    TFR_this_subject=cat(1,TFR_this_subject,Relch);

                end
            end
        end    
        
        if isempty(TFR_this_subject)
            continue
        else
            TFR_subject=nanmean(TFR_this_subject,1);
            TFR_subject=reshape(TFR_subject,1,size(freq.freq,2),size(freq.time,2));
            TFR_this_drug = cat(4,TFR_this_drug, TFR_subject);

              end
    end

    TFRhann_drug           = TFR_pow;
    TFRhann_drug.powspctrm = TFR_this_drug;

   
end

save([Path_variables,'/Figure3/40Hz_PL_EAC_ITPC.mat'],'TFRhann_drug')

%% plot Figure3_2


load([Path_variables,'/Figure3/40Hz_PL_EAC_ITPC.mat'])

TFRhann_drug40_EAC=TFRhann_drug;
TFRhann_drug40_EAC.powspctrm=(mean(TFRhann_drug40_EAC.powspctrm,4));
TFRhann_drug40_EAC.powspctrm=reshape(TFRhann_drug40_EAC.powspctrm,size(TFRhann_drug.freq,2),size(TFRhann_drug.time,2));

clims = [-0.15 0.15];
    
figure
s1 = subplot(1,3,1);
imagesc(TFRhann_drug40_EAC.time,TFRhann_drug40_EAC.freq,TFRhann_drug40_EAC.powspctrm,clims) 
xline(0,'--k')
xline(2,'--k')
yline(40,'--k')
xlim([-0.25 2.225])
xticks([0 1 2])
caxis([0 0.25])
ylabel('Frequency (Hz)')
yticks([5 40 60])
title('Total Power')
set(gca,'YDir','normal')
box off
pbaspect([1 2 1])
ax=gca;
ax.FontSize = 7;


set(s1,'Units','centimeters', 'position', [2.71514378962627 1.62983333333333 2.0113 4.0217]);

cd(Path_figures)

saveas(gcf,'Figure3_2_A.png')