close all

clearvars -except Behav_path Path_figures Path_variables MEG_path

SubjectID={'S01' 'S02' 'S05' 'S15' 'S19' 'S20' 'S22' 'S24','S21','S25','S23'};
Drug=[3 1 1 2 3 2 3; 3 1 3 2 1 2 0; 2 1 1 3 2 3 0; 1 3 3 2 2 1 0; 2 3 1 2 1 3 0; 1 2 3 1 3 2 0; 3 1 2 2 1 3 0; 3 1 1 2 3 2 0; 2 1 1 2 3 3 0; 3 2 1 3 1 2 0;2 3 2 1 1 3 0];
Stimulus=[40 36.75 36.75 40 40 40 40; 40 36.75 36.75 40 40 40 40; 40 36.75 36.75 36.75 40 40 40; 40 36.75 36.75 40 40 40 40; 40 40 40 40 40 40 40; 40 40 40 40 40 40 40; 40 40 40 40 40 40 40; 40 40 40 40 40 40 40;40 40 40 40 40 40 40;40 40 40 40 40 40 40;40 40 40 40 40 40 40];
Freq=[5:0.25:60];
load([Path_variables,'/Figure2/Significant_Channels_stat.mat'])

Significant_channels=[];
if isfield(stat,'posclusterslabelmat')
    for i=1:size(stat.posclusters,2)
        if stat.posclusters(i).prob<0.05
            Significant_channels=unique([Significant_channels;stat.label(find(stat.posclusterslabelmat==i))]);
        end
    end
end

%% Significant Sensors

for Type=2
    TFR_this_drug=[];

    for sub=1:size(SubjectID,2) 
        Sessions = find(Drug(sub,:)==Type);
        StimSessions = Stimulus(sub,Sessions);

        TFR_this_subject=[];

        if isempty(Sessions)==0    

            for Session=1:size(Sessions,2)
                
                % load sensor data
                if (Type==2 && sub==3 && Session==2) 
                    load([MEG_path,'/sensors/',cell2mat(SubjectID(sub)),'/SES',num2str(Sessions(Session)),'/MNE_data_clean_postICA_',cell2mat(SubjectID(sub)),'SES',num2str(Sessions(Session)),'_clean.mat'])
                else
                    load([MEG_path,'/sensors/',cell2mat(SubjectID(sub)),'/SES',num2str(Sessions(Session)),'/MNE_data_clean_postICA_',cell2mat(SubjectID(sub)),'SES',num2str(Sessions(Session)),'.mat'])
                end

                StimFreq = StimSessions(1,Session);

                % 40 Hz sessions
                if StimFreq==40
                
                    id_ripple=find(data.trialinfo(:,1)==1);

                    cfg              = [];
                    cfg.toilim       = [-0.5 2.5-0.0025];
                    data             = ft_redefinetrial(cfg,data);

                    cfg                 = [];
                    cfg.method          = 'template';
                    cfg.template        = 'CTF275_neighb.mat';
                    neighbours          = ft_prepare_neighbours(cfg, data);
                    
                    cfg                 = [];
                    cfg.method          = 'sincos';
                    cfg.neighbours      = neighbours;
                    cfg.channel         = 'MEG';
                    data_planar         = ft_megplanar(cfg, data);
    
                    cfg              = [];
                    cfg.output       = 'pow';
                    cfg.channel      = Significant_channels;
                    cfg.trials       = id_ripple;
                    cfg.method       = 'mtmconvol';
                    cfg.taper        = 'hanning';
                    cfg.foi          = Freq;      
                    cfg.pad          = 4;
                    cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   
                    cfg.toi          = -0.5:0.025:(2.5-0.0025); 
                    TFR_this_session = ft_freqanalysis(cfg, data_planar);
    
                    % relative power change  to baseline 
                    cfg              = [];
                    cfg.baseline     = [-0.5 0];
                    cfg.baselinetype = 'relchange';
                    Relch            = ft_freqbaseline(cfg,TFR_this_session);
    
                    TFR_this_subject=cat(4,TFR_this_subject,Relch.powspctrm);
                end
            end
        end
        
        % concatenate across subjects
        TFR_this_subject=mean(TFR_this_subject,4);
        TFR_this_drug = cat(4,TFR_this_drug,TFR_this_subject);
    
    end
    
    
    TFR_drug_B=TFR_this_session;
    TFR_drug_B.powspctrm=TFR_this_drug;
end

save([Path_variables,'/Figure2/40Hz_PL_Sensor_sig.mat'],'TFR_drug_B')

%% All Sensors

for Type=2

    TFR_this_drug=[];

    for sub=1:size(SubjectID,2) 
        Sessions = find(Drug(sub,:)==Type);
        StimSessions = Stimulus(sub,Sessions);

        TFR_this_subject=[];

        if isempty(Sessions)==0    

            for Session=1:size(Sessions,2)

                if (Type==2 && sub==3 && Session==2) 
                    load([MEG_path,'/sensors/',cell2mat(SubjectID(sub)),'/SES',num2str(Sessions(Session)),'/MNE_data_clean_postICA_',cell2mat(SubjectID(sub)),'SES',num2str(Sessions(Session)),'_clean.mat'])
                else
                    load([MEG_path,'/sensors/',cell2mat(SubjectID(sub)),'/SES',num2str(Sessions(Session)),'/MNE_data_clean_postICA_',cell2mat(SubjectID(sub)),'SES',num2str(Sessions(Session)),'.mat'])
                end

                StimFreq = StimSessions(1,Session);

                if StimFreq==40
                
                    id_ripple=find(data.trialinfo(:,1)==1);

                    cfg              = [];
                    cfg.toilim       = [-0.5 2.5-0.0025];
                    data             = ft_redefinetrial(cfg,data);

                    cfg                 = [];
                    cfg.method          = 'template';
                    cfg.template        = 'CTF275_neighb.mat';
                    neighbours          = ft_prepare_neighbours(cfg, data);
                    
                    cfg                 = [];
                    cfg.method          = 'sincos';
                    cfg.neighbours      = neighbours;
                    data_planar         = ft_megplanar(cfg, data);
    
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
                    TFR_this_session = ft_freqanalysis(cfg, data_planar);
    
                    cfg                 = [];
                    TFR_this_session    = ft_combineplanar(cfg, TFR_this_session);
    
                    cfg              = [];
                    cfg.baseline     = [-0.5 0];
                    cfg.baselinetype = 'relchange';
                    Relch            = ft_freqbaseline(cfg,TFR_this_session);
    
                    TFR_this_subject=cat(4,TFR_this_subject,Relch.powspctrm);
                end
            end
        end
        TFR_this_subject=mean(TFR_this_subject,4);
        TFR_this_drug = cat(4,TFR_this_drug,TFR_this_subject);
    end
    TFR_drug_B=TFR_this_session;
    TFR_drug_B.powspctrm=TFR_this_drug;
end

save([Path_variables,'/Figure2/40Hz_PL_Sensor_all.mat'],'TFR_drug_B')