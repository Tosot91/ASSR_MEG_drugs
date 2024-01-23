clearvars -except Behav_path Path_figures Path_variables

close all

% caculate TFR separately  for each drug

SubjectID={'S01' 'S02' 'S04' 'S05' 'S06' 'S07' 'S08' 'S10' 'S11' 'S12' 'S15' 'S17' 'S18' 'S19' 'S20' 'S22' 'S24','S21','S25','S23'};
Drug=[3 1 1 2 3 2 3; 3 1 3 2 1 2 0; 1 3 2 1 3 2 0; 2 1 1 3 2 3 0; 2 1 3 1 2 3 0; 1 3 2 3 1 2 0; 3 1 1 2 2 3 0; 1 2 2 3 1 3 1; 1 2 2 3 3 1 1; 3 1 2 1 3 2 0; 1 3 3 2 2 1 0; 2 1 3 2 3 1 2; 1 2 3 2 1 3 1; 2 3 1 2 1 3 0; 1 2 3 1 3 2 0; 3 1 2 2 1 3 0; 3 1 1 2 3 2 0; 2 1 1 2 3 3 0; 3 2 1 3 1 2 0;2 3 2 1 1 3 0];
Stimulus=[40 36.75 36.75 40 40 40 40; 40 36.75 36.75 40 40 40 40; 40 36.75 36.75 36.75 40 40 40; 40 36.75 36.75 36.75 40 40 40; 36.75 40 40 40 40 40 40; 40 36.75 36.75 36.75 40 40 40; 40 36.75 36.75 36.75 40 40 40; 40 36.75 36.75 40 40 40 40; 40 36.75 36.75 36.75 40 40 40; 40 36.75 36.75 36.75 40 40 40; 40 36.75 36.75 40 40 40 40; 40 36.75 36.75 36.75 40 40 40; 40 36.75 36.75 36.75 40 40 40; 40 40 40 40 40 40 40; 40 40 40 40 40 40 40; 40 40 40 40 40 40 40; ...
    40 40 40 40 40 40 40;40 40 40 40 40 40 40;40 40 40 40 40 40 40;40 40 40 40 40 40 40];
Freq=[29:0.25:48];

% drug type: 1 LO, 2 PL, 3 ME ;
for Type=1:3
    TFR_this_drug=[];
    TFR_this_drug_noBase=[];
    TFR_this_drug_phase=[];
    TFR_this_drug_phase_noBase=[];

    for sub=1:size(SubjectID,2) 
        Sessions = find(Drug(sub,:)==Type);
        StimSessions = Stimulus(sub,Sessions);
    
        TFR_this_subject=[];
        TFR_this_subject_noBase=[];
        TFR_this_subject_phase=[];
        TFR_this_subject_phase_noBase=[];
    
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
                
    
                data.trial=datas_all;
                data.label=data.label(1:size(datas_all,2),:);
    
                StimFreq = StimSessions(1,Session);
    
                id_ripple=find(data.trialinfo(:,1)==1);
    
                % time-series for TFR of power change
                cfg              = [];
                cfg.toilim       = [-0.5 2.5-0.0025];
                data             = ft_redefinetrial(cfg,data);
                
                % ERF for TFR of phase-locked analyses
   
                cfg=[];
                cfg.trials = id_ripple;
                ERF = ft_timelockanalysis(cfg,data);
                 % TFR  phase-locked 
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
    
                % get -7/+7  Hz around stiulus frequency
               
                [val,idx]=min(abs(TFR_phaselocked.freq-StimFreq));
                TFR_phaselocked.powspctrm=TFR_phaselocked.powspctrm(:,:,idx-28:idx+28,:);
                TFR_phaselocked.freq=-7:0.25:7;
                
                TFR_phaselocked.powspctrm=mean(TFR_phaselocked.powspctrm,2);
                TFR_phaselocked.label={'EAC'};
                
                % TFR  total power
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
                    
                % get -7/+7  Hz around stiulus frequency
                [val,idx]=min(abs(TFR_this_session.freq-StimFreq));
    
                TFR_this_session.powspctrm=TFR_this_session.powspctrm(:,:,idx-28:idx+28,:);
                TFR_this_session.freq=-7:0.25:7;
                TFR_this_session.powspctrm=mean(TFR_this_session.powspctrm, 2);
                TFR_this_session.label={'EAC'};
                
                % matrix of pre-stimulus baseline 
                Baseline_this_session=nanmean(nanmean(nanmean(TFR_this_session.powspctrm(:,:,:,1:21),4),2),1);
                Baseline_this_session=reshape(Baseline_this_session,1,size(TFR_this_session.freq,2));
                Base_all=repmat(Baseline_this_session,[size(id_ripple,1),1,size(TFR_this_session.time,2)]);
                Base_all=reshape(Base_all,size(id_ripple,1),1,size(TFR_this_session.freq,2),size(TFR_this_session.time,2));
    
                Baseline_phaselocked=nanmean(nanmean(nanmean(TFR_phaselocked.powspctrm(:,:,:,1:21),4),2),1);
                Baseline_phaselocked=reshape(Baseline_phaselocked,1,size(TFR_phaselocked.freq,2));
                Base_all_phaselocked=repmat(Baseline_phaselocked,[1,1,size(TFR_phaselocked.time,2)]);
                Base_all_phaselocked=reshape(Base_all_phaselocked,1,1,size(TFR_phaselocked.freq,2),size(TFR_phaselocked.time,2));
                Base_all_phase=repmat(Baseline_this_session,[1,1,size(TFR_phaselocked.time,2)]);
                Base_all_phase=reshape(Base_all_phase,1,1,size(TFR_phaselocked.freq,2),size(TFR_phaselocked.time,2));
                
                
                % get relative change in power 
                Relch = (TFR_this_session.powspctrm-Base_all)./Base_all;
                Relch_phase = (TFR_phaselocked.powspctrm-Base_all_phaselocked)./Base_all_phase;
                
                % concatenate across sessions  for this  subject
                TFR_this_subject=cat(1,TFR_this_subject,Relch);
                TFR_this_subject_noBase=cat(1,TFR_this_subject_noBase,TFR_this_session.powspctrm);
                TFR_this_subject_phase=cat(1,TFR_this_subject_phase,Relch_phase);
                TFR_this_subject_phase_noBase=cat(1,TFR_this_subject_phase_noBase,TFR_phaselocked.powspctrm);
            end
            
            % average across sessions anc cocatenate across subjects
            TFR_this_subject=mean(TFR_this_subject,1);
            TFR_this_drug = cat(1,TFR_this_drug,TFR_this_subject);

            TFR_this_subject_noBase=mean(TFR_this_subject_noBase,1);
            TFR_this_drug_noBase = cat(1,TFR_this_drug_noBase,TFR_this_subject_noBase);
           
            TFR_this_subject_phase=mean(TFR_this_subject_phase,1);
            TFR_this_drug_phase = cat(1,TFR_this_drug_phase,TFR_this_subject_phase);

            TFR_this_subject_phase_noBase=mean(TFR_this_subject_phase_noBase,1);
            TFR_this_drug_phase_noBase = cat(1,TFR_this_drug_phase_noBase,TFR_this_subject_phase_noBase);
        end
    end
 
    
    % save individual Drugs results as separate structures
    if Type==1
        TFR_drug_A=TFR_this_session;
        TFR_drug_A.dimord='subj_chan_freq_time';
        TFR_drug_A.powspctrm=TFR_this_drug;

        TFR_drug_A_noBase=TFR_this_session;
        TFR_drug_A_noBase.dimord='subj_chan_freq_time';
        TFR_drug_A_noBase.powspctrm=TFR_this_drug_noBase;

        TFR_drug_A_phase=TFR_this_session;
        TFR_drug_A_phase.dimord='subj_chan_freq_time';
        TFR_drug_A_phase.powspctrm=TFR_this_drug_phase;

        TFR_drug_A_phase_noBase=TFR_this_session;
        TFR_drug_A_phase_noBase.dimord='subj_chan_freq_time';
        TFR_drug_A_phase_noBase.powspctrm=TFR_this_drug_phase_noBase;

    elseif Type==2
        TFR_drug_B=TFR_this_session;
        TFR_drug_B.dimord='subj_chan_freq_time';
        TFR_drug_B.powspctrm=TFR_this_drug;

        TFR_drug_B_noBase=TFR_this_session;
        TFR_drug_B_noBase.dimord='subj_chan_freq_time';
        TFR_drug_B_noBase.powspctrm=TFR_this_drug_noBase;

        TFR_drug_B_phase=TFR_this_session;
        TFR_drug_B_phase.dimord='subj_chan_freq_time';
        TFR_drug_B_phase.powspctrm=TFR_this_drug_phase;

        TFR_drug_B_phase_noBase=TFR_this_session;
        TFR_drug_B_phase_noBase.dimord='subj_chan_freq_time';
        TFR_drug_B_phase_noBase.powspctrm=TFR_this_drug_phase_noBase;

    elseif Type==3
        TFR_drug_C=TFR_this_session;
        TFR_drug_C.dimord='subj_chan_freq_time';
        TFR_drug_C.powspctrm=TFR_this_drug;

        TFR_drug_C_noBase=TFR_this_session;
        TFR_drug_C_noBase.dimord='subj_chan_freq_time';
        TFR_drug_C_noBase.powspctrm=TFR_this_drug_noBase;

        TFR_drug_C_phase=TFR_this_session;
        TFR_drug_C_phase.dimord='subj_chan_freq_time';
        TFR_drug_C_phase.powspctrm=TFR_this_drug_phase;

        TFR_drug_C_phase_noBase=TFR_this_session;
        TFR_drug_C_phase_noBase.dimord='subj_chan_freq_time';
        TFR_drug_C_phase_noBase.powspctrm=TFR_this_drug_phase_noBase;
    end
end
%%
save([Path_variables,'/Figure3/GrandAvg_drug_A_EAC.mat'],'TFR_drug_A')

save([Path_variables,'/Figure3/GrandAvg_drug_A_EAC_noBase.mat'],'TFR_drug_A_noBase')

save([Path_variables,'/Figure3/GrandAvg_drug_A_EAC_phase.mat'],'TFR_drug_A_phase')

save([Path_variables,'/Figure3/GrandAvg_drug_A_EAC_phase_noBase.mat'],'TFR_drug_A_phase_noBase')

save([Path_variables,'/Figure3/GrandAvg_drug_B_EAC.mat'],'TFR_drug_B')

save([Path_variables,'/Figure3/GrandAvg_drug_B_EAC_noBase.mat'],'TFR_drug_B_noBase')

save([Path_variables,'/Figure3/GrandAvg_drug_B_EAC_phase.mat'],'TFR_drug_B_phase')

save([Path_variables,'/Figure3/GrandAvg_drug_B_EAC_phase_noBase.mat'],'TFR_drug_B_phase_noBase')

save([Path_variables,'/Figure3/GrandAvg_drug_C_EAC.mat'],'TFR_drug_C')

save([Path_variables,'/Figure3/GrandAvg_drug_C_EAC_noBase.mat'],'TFR_drug_C_noBase')

save([Path_variables,'/Figure3/GrandAvg_drug_C_EAC_phase.mat'],'TFR_drug_C_phase')

save([Path_variables,'/Figure3/GrandAvg_drug_C_EAC_phase_noBase.mat'],'TFR_drug_C_phase_noBase')