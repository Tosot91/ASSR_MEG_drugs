close all

clearvars -except Behav_path Path_figures Path_variables Code_path MEG_path



SubjectID={'S01' 'S02' 'S04' 'S05' 'S06' 'S07' 'S08' 'S10' 'S11' 'S12' 'S15' 'S17' 'S18' 'S19' 'S20' 'S22' 'S24','S21','S25','S23'};
Drug=[3 1 1 2 3 2 3; 3 1 3 2 1 2 0; 1 3 2 1 3 2 0; 2 1 1 3 2 3 0; 2 1 3 1 2 3 0; 1 3 2 3 1 2 0; 3 1 1 2 2 3 0; 1 2 2 3 1 3 1; 1 2 2 3 3 1 1; 3 1 2 1 3 2 0; 1 3 3 2 2 1 0; 2 1 3 2 3 1 2; 1 2 3 2 1 3 1; 2 3 1 2 1 3 0; 1 2 3 1 3 2 0; 3 1 2 2 1 3 0; 3 1 1 2 3 2 0; 2 1 1 2 3 3 0; 3 2 1 3 1 2 0;2 3 2 1 1 3 0];
Stimulus=[40 36.75 36.75 40 40 40 40; 40 36.75 36.75 40 40 40 40; 40 36.75 36.75 36.75 40 40 40; 40 36.75 36.75 36.75 40 40 40; 36.75 40 40 40 40 40 40; 40 36.75 36.75 36.75 40 40 40; 40 36.75 36.75 36.75 40 40 40; 40 36.75 36.75 40 40 40 40; 40 36.75 36.75 36.75 40 40 40; 40 36.75 36.75 36.75 40 40 40; 40 36.75 36.75 40 40 40 40; 40 36.75 36.75 36.75 40 40 40; 40 36.75 36.75 36.75 40 40 40; 40 40 40 40 40 40 40; 40 40 40 40 40 40 40; 40 40 40 40 40 40 40; ...
    40 40 40 40 40 40 40;40 40 40 40 40 40 40;40 40 40 40 40 40 40;40 40 40 40 40 40 40];
Freq=[29:0.25:48];
% calculate ITPC for all drugs

for Type=1:3
    
    ITC_this_drug=[];
    ITC_this_drug_noBase=[];

    for sub=1:size(SubjectID,2) 
        
        Sessions = find(Drug(sub,:)==Type);
        StimSessions = Stimulus(sub,Sessions);
    
        ITC_this_subject=[];
        ITC_this_subject_noBase=[];
    
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
    
                cfg              = [];
                cfg.toilim       = [-0.5 2.5-0.0025];
                data             = ft_redefinetrial(cfg,data);
% 
                if sub==1
                    cfg              = [];
                    cfg.output       = 'pow';
                    cfg.channel      = 'MEG';
                    cfg.trials       = id_ripple;
                    cfg.method       = 'mtmconvol';
                    cfg.taper        = 'hanning';
                    cfg.foi          = Freq;     
                    cfg.pad          = 4;
                    cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   
                    cfg.toi          = -0.5:0.025:(2.5-0.0025);                 
                    TFR_this_session = ft_freqanalysis(cfg, data);
                        
                    [val,idx]=min(abs(TFR_this_session.freq-StimFreq));
        
                    TFR_this_session.powspctrm=TFR_this_session.powspctrm(:,idx-28:idx+28,:);
                    TFR_this_session.freq=-7:0.25:7;
                    TFR_this_session.powspctrm=mean(TFR_this_session.powspctrm, 2);
                    TFR_this_session.label={'EAC'};
                end

                
                % calculate ITPC  at +/- 7 Hz from stimulus frequency
                
                    cfg              = [];
                    cfg.output       = 'fourier';
                    cfg.channel      = 'MEG';
                    cfg.trials       = id_ripple;
                    cfg.method       = 'mtmconvol';
                    cfg.taper        = 'hanning';
                    cfg.foi          = Freq;     
                    cfg.pad          = 4;
                    cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   
                    cfg.toi          = -0.5:0.025:(2.5-0.0025);                 
                    freq = ft_freqanalysis(cfg, data);
                
                [val,idx]=min(abs(freq.freq-StimFreq));
    
                freq.fourierspctrm=freq.fourierspctrm(:,:,idx-28:idx+28,:);
                freq.freq=-7:0.25:7;
                
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

    
                ITC_this_subject        = cat(4,ITC_this_subject,itc.itpc); 
            end
        end
        ITC_this_subject=mean(ITC_this_subject,4);
        ITC_this_drug = cat(4,ITC_this_drug,ITC_this_subject);

    end
     
    % save different files for different drug types

    if Type==1
        ITC_drug_A=TFR_this_session;
         ITC_drug_A.powspctrm=ITC_this_drug;


save([Path_variables,'/Figure3/GrandAvg_drug_A_EAC_ITPC.mat'],'ITC_drug_A')

    elseif Type==2
         ITC_drug_B=TFR_this_session;
         ITC_drug_B.powspctrm=ITC_this_drug;

save([Path_variables,'/Figure3/GrandAvg_drug_B_EAC_ITPC.mat'],'ITC_drug_B')
    elseif Type==3
         ITC_drug_C=TFR_this_session;
         ITC_drug_C.powspctrm=ITC_this_drug;

save([Path_variables,'/Figure3/GrandAvg_drug_C_EAC_ITPC.mat'],'ITC_drug_C')
    end
end