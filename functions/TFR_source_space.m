clearvars -except Behav_path Path_figures Path_variables MEG_path

close all

SubjectID={'S01' 'S02' 'S05' 'S15' 'S19' 'S20' 'S22' 'S24','S21','S25','S23'};
Drug=[3 1 1 2 3 2 3; 3 1 3 2 1 2 0; 2 1 1 3 2 3 0; 1 3 3 2 2 1 0; 2 3 1 2 1 3 0; 1 2 3 1 3 2 0; 3 1 2 2 1 3 0; 3 1 1 2 3 2 0; 2 1 1 2 3 3 0; 3 2 1 3 1 2 0;2 3 2 1 1 3 0];
Stimulus=[40 36.75 36.75 40 40 40 40; 40 36.75 36.75 40 40 40 40; 40 36.75 36.75 36.75 40 40 40; 40 36.75 36.75 40 40 40 40; 40 40 40 40 40 40 40; 40 40 40 40 40 40 40; 40 40 40 40 40 40 40; 40 40 40 40 40 40 40;40 40 40 40 40 40 40;40 40 40 40 40 40 40;40 40 40 40 40 40 40];
Freq=[5:0.25:60];

%% 40Hz EAC

% TFR for  40Hz AM placebo session,  at source level (EAC)
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

                if StimFreq==40

                    id_ripple = find(data.trialinfo(:,1)==1);

                    % data for TFR power change
                    cfg              = [];
                    cfg.toilim       = [-0.5 2.5-0.0025];
                    data             = ft_redefinetrial(cfg,data);

                    % ERF averaged across ripple trials
                    cfg              = [];
                    cfg.trials       = id_ripple;
                    ERF              = ft_timelockanalysis(cfg,data);

                    %TFR power change
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
                    TFR_this_session = ft_freqanalysis(cfg, data);

                    % TFR phase-locked
                    cfg              = [];
                    cfg.output       = 'pow';
                    cfg.channel      = 'MEG';
                    cfg.keeptrials   = 'yes';
                    cfg.method       = 'mtmconvol';
                    cfg.taper        = 'hanning';
                    cfg.foi          = Freq;     
                    cfg.pad          = 4;
                    cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   
                    cfg.toi          = -0.5:0.025:(2.5-0.0025);             
                    TFR_phaselocked  = ft_freqanalysis(cfg, ERF);
                   
                    TFR_this_session.powspctrm=mean(TFR_this_session.powspctrm, 2);
                    TFR_this_session.label={'EAC'};
                    
                    % baseline -0.5 to 0
                    Baseline_this_session=nanmean(nanmean(nanmean(TFR_this_session.powspctrm(:,:,:,1:21),4),2),1);
                    Baseline_this_session=reshape(Baseline_this_session,1,size(TFR_this_session.freq,2));
                   
                    Base_all=repmat(Baseline_this_session,[size(id_ripple,1),1,size(TFR_this_session.time,2)]);
                    Base_all=reshape(Base_all,size(id_ripple,1),1,size(TFR_this_session.freq,2),size(TFR_this_session.time,2));
                    
                    % relative change in power 
                    Relch = (TFR_this_session.powspctrm-Base_all)./Base_all;
                    TFR_this_subject=cat(1,TFR_this_subject,Relch);
                    
                    
                    TFR_phaselocked.powspctrm=mean(TFR_phaselocked.powspctrm, 2);
                    TFR_phaselocked.label={'EAC'};
                    
                    % same baseline for phase-locked analyses
                    Baseline_phaselocked=nanmean(nanmean(nanmean(TFR_phaselocked.powspctrm(:,:,:,1:21),4),2),1);
                    Baseline_phaselocked=reshape(Baseline_phaselocked,1,size(TFR_this_session.freq,2));
                    Base_all_phaselocked=repmat(Baseline_phaselocked,[1,1,size(TFR_this_session.time,2)]);
                    Base_all_phaselocked=reshape(Base_all_phaselocked,1,1,size(TFR_this_session.freq,2),size(TFR_this_session.time,2));
                    Base_all_phase=repmat(Baseline_this_session,[1,1,size(TFR_this_session.time,2)]);
                    Base_all_phase=reshape(Base_all_phase,1,1,size(TFR_this_session.freq,2),size(TFR_this_session.time,2));
                    
                    % relch phase-locked power
                    Relch_phase = (TFR_phaselocked.powspctrm-Base_all_phaselocked)./Base_all_phase;
                    TFR_this_subject_phase=cat(1,TFR_this_subject_phase,Relch_phase);

                end
            end
        end    
        
        if isempty(TFR_this_subject)
            continue
        else
            
            %concatenate across subejcts
            TFR_subject=nanmean(TFR_this_subject,1);
            TFR_subject=reshape(TFR_subject,1,size(TFR_this_session.freq,2),size(TFR_this_session.time,2));
            TFR_this_drug = cat(4,TFR_this_drug, TFR_subject);

            TFR_subject_phaselocked=nanmean(TFR_this_subject_phase,1);
            TFR_subject_phaselocked=reshape(TFR_subject_phaselocked,1,size(TFR_this_session.freq,2),size(TFR_this_session.time,2));
            TFR_this_drug_phase = cat(4,TFR_this_drug_phase, TFR_subject_phaselocked);
        end
    end

    TFRhann_drug           = TFR_this_session;
    TFRhann_drug.powspctrm = TFR_this_drug;

    TFRhann_drug_phase           = TFR_this_session;
    TFRhann_drug_phase.powspctrm = TFR_this_drug_phase;

end

save([Path_variables,'/Figure2/40Hz_PL_EAC.mat'],'TFRhann_drug')

save([Path_variables,'/Figure2/40Hz_PL_phase_EAC.mat'],'TFRhann_drug_phase')

%% same as before but for dlPFC

for Type=2
  
    TFR_this_drug=[];
     
    for sub=1:size(SubjectID,2)
        Sessions = find(Drug(sub,:)==Type);
        StimSessions = Stimulus(sub,Sessions);

        TFR_this_subject=[];

        if isempty(Sessions)==0
            
            for Session=1:size(Sessions,2)

                if Session==3
                    load([MEG_path,'/source/',cell2mat(SubjectID(sub)),'/SES11/9-46d_ERF_l_long.mat'])
                    datas1=comb_dict{1,1};
                    load([MEG_path,'/source/',cell2mat(SubjectID(sub)),'/SES11/p9-46v_ERF_l_long.mat'])
                    datas2=comb_dict{1,1};
                    load([MEG_path,'/source/',cell2mat(SubjectID(sub)),'/SES11/a9-46v_ERF_l_long.mat'])
                    datas3=comb_dict{1,1};
                    load([MEG_path,'/source/',cell2mat(SubjectID(sub)),'/SES11/46_ERF_l_long.mat'])
                    datas4=comb_dict{1,1};
                    load([MEG_path,'/source/',cell2mat(SubjectID(sub)),'/SES11/9-46d_ERF_r_long.mat'])
                    datas5=comb_dict{1,1};
                    load([MEG_path,'/source/',cell2mat(SubjectID(sub)),'/SES11/p9-46v_ERF_r_long.mat'])
                    datas6=comb_dict{1,1};
                    load([MEG_path,'/source/',cell2mat(SubjectID(sub)),'/SES11/a9-46v_ERF_r_long.mat'])
                    datas7=comb_dict{1,1};
                    load([MEG_path,'/source/',cell2mat(SubjectID(sub)),'/SES11/46_ERF_r_long.mat'])
                    datas8=comb_dict{1,1};
                else
                    load([MEG_path,'/source/',cell2mat(SubjectID(sub)),'/SES',num2str(Sessions(Session)),'/9-46d_ERF_l_long.mat'])
                    datas1=comb_dict{1,1};
                    load([MEG_path,'/source/',cell2mat(SubjectID(sub)),'/SES',num2str(Sessions(Session)),'/p9-46v_ERF_l_long.mat'])
                    datas2=comb_dict{1,1};
                    load([MEG_path,'/source/',cell2mat(SubjectID(sub)),'/SES',num2str(Sessions(Session)),'/a9-46v_ERF_l_long.mat'])
                    datas3=comb_dict{1,1};
                    load([MEG_path,'/source/',cell2mat(SubjectID(sub)),'/SES',num2str(Sessions(Session)),'/46_ERF_l_long.mat'])
                    datas4=comb_dict{1,1};
                    load([MEG_path,'/source/',cell2mat(SubjectID(sub)),'/SES',num2str(Sessions(Session)),'/9-46d_ERF_r_long.mat'])
                    datas5=comb_dict{1,1};
                    load([MEG_path,'/source/',cell2mat(SubjectID(sub)),'/SES',num2str(Sessions(Session)),'/p9-46v_ERF_r_long.mat'])
                    datas6=comb_dict{1,1};
                    load([MEG_path,'/source/',cell2mat(SubjectID(sub)),'/SES',num2str(Sessions(Session)),'/a9-46v_ERF_r_long.mat'])
                    datas7=comb_dict{1,1};
                    load([MEG_path,'/source/',cell2mat(SubjectID(sub)),'/SES',num2str(Sessions(Session)),'/46_ERF_r_long.mat'])
                    datas8=comb_dict{1,1};
                end

                datas_all = cat(2,datas1.erfdata,datas2.erfdata,datas3.erfdata,datas4.erfdata,datas5.erfdata,datas6.erfdata,datas7.erfdata,datas8.erfdata);

                if Session==3
                    load([MEG_path,'/sensors/',cell2mat(SubjectID(sub)),'/SES11/MNE_data_clean_postICA_',cell2mat(SubjectID(sub)),'SES11.mat'])
                else
                    load([MEG_path,'/sensors/',cell2mat(SubjectID(sub)),'/SES',num2str(Sessions(Session)),'/MNE_data_clean_postICA_',cell2mat(SubjectID(sub)),'SES',num2str(Sessions(Session)),'.mat'])
                end
                
                data.trial = datas_all;
                data.label = data.label(1:size(datas_all,2),:);
    
                StimFreq = StimSessions(1,Session);

                if StimFreq==40

                    id_ripple = find(data.trialinfo(:,1)==1);

                    cfg              = [];
                    cfg.toilim       = [-0.5 2.5-0.0025];
                    data             = ft_redefinetrial(cfg,data);

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
                    TFR_this_session = ft_freqanalysis(cfg, data);

                    TFR_this_session.powspctrm=mean(TFR_this_session.powspctrm, 2);
                    TFR_this_session.label={'EAC'};
                    Baseline_this_session=nanmean(nanmean(nanmean(TFR_this_session.powspctrm(:,:,:,1:21),4),2),1);
                    Baseline_this_session=reshape(Baseline_this_session,1,size(TFR_this_session.freq,2));
                    Base_all=repmat(Baseline_this_session,[size(id_ripple,1),1,size(TFR_this_session.time,2)]);
                    Base_all=reshape(Base_all,size(id_ripple,1),1,size(TFR_this_session.freq,2),size(TFR_this_session.time,2));
                    Relch = (TFR_this_session.powspctrm-Base_all)./Base_all;
                    TFR_this_subject=cat(1,TFR_this_subject,Relch);
                end
            end
        end    
        
        if isempty(TFR_this_subject)
            continue
        else
            TFR_subject=nanmean(TFR_this_subject,1);
            TFR_subject=reshape(TFR_subject,1,size(TFR_this_session.freq,2),size(TFR_this_session.time,2));
            TFR_this_drug = cat(4,TFR_this_drug, TFR_subject);
        end

    end
    TFRhann_drug           = TFR_this_session;
    TFRhann_drug.powspctrm = TFR_this_drug;
end

save([Path_variables,'/Figure2/40Hz_PL_dlPFC.mat'],'TFRhann_drug')

%% same as before but for 37-hz AM stimulus session

SubjectID={'S10' 'S11' 'S18'};
Drug=[1 2 2 3 1 3 1; 1 2 2 3 3 1 1; 1 2 3 2 1 3 1];
Stimulus=[40 36.75 36.75 40 40 40 40; 40 36.75 36.75 36.75 40 40 40; 40 36.75 36.75 36.75 40 40 40];
Freq=[5:0.25:60];

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

                if StimFreq==36.75

                    id_ripple = find(data.trialinfo(:,1)==1);

                    cfg              = [];
                    cfg.toilim       = [-0.5 2.5-0.0025];
                    data             = ft_redefinetrial(cfg,data);

                    cfg=[];
                    cfg.trials = id_ripple;
                    ERF = ft_timelockanalysis(cfg,data);

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
                    TFR_this_session = ft_freqanalysis(cfg, data);

                    cfg              = [];
                    cfg.output       = 'pow';
                    cfg.channel      = 'MEG';
                    cfg.keeptrials   = 'yes';
                    cfg.method       = 'mtmconvol';
                    cfg.taper        = 'hanning';
                    cfg.foi          = Freq;     
                    cfg.pad          = 4;
                    cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   
                    cfg.toi          = -0.5:0.025:(2.5-0.0025);             
                    TFR_phaselocked  = ft_freqanalysis(cfg, ERF);

                    TFR_this_session.powspctrm=mean(TFR_this_session.powspctrm, 2);
                    TFR_this_session.label={'EAC'};
                    Baseline_this_session=nanmean(nanmean(nanmean(TFR_this_session.powspctrm(:,:,:,1:21),4),2),1);
                    Baseline_this_session=reshape(Baseline_this_session,1,size(TFR_this_session.freq,2));
                    Base_all=repmat(Baseline_this_session,[size(id_ripple,1),1,size(TFR_this_session.time,2)]);
                    Base_all=reshape(Base_all,size(id_ripple,1),1,size(TFR_this_session.freq,2),size(TFR_this_session.time,2));
                    Relch = (TFR_this_session.powspctrm-Base_all)./Base_all;
                    TFR_this_subject=cat(1,TFR_this_subject,Relch);

                    TFR_phaselocked.powspctrm=mean(TFR_phaselocked.powspctrm, 2);
                    TFR_phaselocked.label={'EAC'};
                    Baseline_phaselocked=nanmean(nanmean(nanmean(TFR_phaselocked.powspctrm(:,:,:,1:21),4),2),1);
                    Baseline_phaselocked=reshape(Baseline_phaselocked,1,size(TFR_this_session.freq,2));
                    Base_all_phaselocked=repmat(Baseline_phaselocked,[1,1,size(TFR_this_session.time,2)]);
                    Base_all_phaselocked=reshape(Base_all_phaselocked,1,1,size(TFR_this_session.freq,2),size(TFR_this_session.time,2));
                    Base_all_phase=repmat(Baseline_this_session,[1,1,size(TFR_this_session.time,2)]);
                    Base_all_phase=reshape(Base_all_phase,1,1,size(TFR_this_session.freq,2),size(TFR_this_session.time,2));
                    Relch_phase = (TFR_phaselocked.powspctrm-Base_all_phaselocked)./Base_all_phase;
                    TFR_this_subject_phase=cat(1,TFR_this_subject_phase,Relch_phase);

                end
            end
        end    
        
        if isempty(TFR_this_subject)
            continue
        else
            TFR_subject=nanmean(TFR_this_subject,1);
            TFR_subject=reshape(TFR_subject,1,size(TFR_this_session.freq,2),size(TFR_this_session.time,2));
            TFR_this_drug = cat(4,TFR_this_drug, TFR_subject);

            TFR_subject_phaselocked=nanmean(TFR_this_subject_phase,1);
            TFR_subject_phaselocked=reshape(TFR_subject_phaselocked,1,size(TFR_this_session.freq,2),size(TFR_this_session.time,2));
            TFR_this_drug_phase = cat(4,TFR_this_drug_phase, TFR_subject_phaselocked);
        end

    end

    TFRhann_drug           = TFR_this_session;
    TFRhann_drug.powspctrm = TFR_this_drug;

    TFRhann_drug_phase           = TFR_this_session;
    TFRhann_drug_phase.powspctrm = TFR_this_drug_phase;

end

save([Path_variables,'/Figure2/37Hz_PL_EAC.mat'],'TFRhann_drug')

save([Path_variables,'/Figure2/37Hz_PL_phase_EAC.mat'],'TFRhann_drug_phase')