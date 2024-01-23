close all

clearvars -except Behav_path Path_figures Path_variables MEG_path

% find sensors at which ASSR is significant for Placebo sessions 

SubjectID={'S01' 'S02' 'S05' 'S15' 'S19' 'S20' 'S22' 'S24','S21','S25','S23'};
Drug=[3 1 1 2 3 2 3; 3 1 3 2 1 2 0; 2 1 1 3 2 3 0; 1 3 3 2 2 1 0; 2 3 1 2 1 3 0; 1 2 3 1 3 2 0; 3 1 2 2 1 3 0; 3 1 1 2 3 2 0; 2 1 1 2 3 3 0; 3 2 1 3 1 2 0;2 3 2 1 1 3 0];
Stimulus=[40 36.75 36.75 40 40 40 40; 40 36.75 36.75 40 40 40 40; 40 36.75 36.75 36.75 40 40 40; 40 36.75 36.75 40 40 40 40; 40 40 40 40 40 40 40; 40 40 40 40 40 40 40; 40 40 40 40 40 40 40; 40 40 40 40 40 40 40;40 40 40 40 40 40 40;40 40 40 40 40 40 40;40 40 40 40 40 40 40];
Freq=[29:0.25:48];

% drug_type:2 is placebo
for Type=2
    
    for sub=1:size(SubjectID,2)
        
        % get sessions  from this subject and this drug type
        Sessions=find(Drug(sub,:)==Type);
        StimSessions = Stimulus(sub,Sessions);
        
        if isempty(Sessions)==0
            Significant_channels=[];
            
            for Session=1:size(Sessions,2)
                
                
                if (Type==2 && sub==3 && Session==2) 
                    load([MEG_path,'/sensors/',cell2mat(SubjectID(sub)),'/SES',num2str(Sessions(Session)),'/MNE_data_clean_postICA_',cell2mat(SubjectID(sub)),'SES',num2str(Sessions(Session)),'_clean.mat'])
                else
                    load([MEG_path,'/sensors/',cell2mat(SubjectID(sub)),'/SES',num2str(Sessions(Session)),'/MNE_data_clean_postICA_',cell2mat(SubjectID(sub)),'SES',num2str(Sessions(Session)),'.mat'])
                end
                
                
                StimFreq = StimSessions(1,Session);
    

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
                cfg.channel      = 'MEG';
                cfg.trials       = id_ripple;
                cfg.method       = 'mtmconvol';
                cfg.taper        = 'hanning';
                cfg.foi          = Freq;      
                cfg.pad          = 4;
                cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   
                cfg.toi          = -0.5:0.025:(2.5-0.0025); 
                cfg.keeptrials   = 'yes';
                TFR_this_session = ft_freqanalysis(cfg, data_planar);
    
                % concatenate acrosssubjects
                if sub==1 && Session==1
                    data_all=TFR_this_session;
                else
                    data_all.powspctrm=cat(1,data_all.powspctrm,TFR_this_session.powspctrm);
                    data_all.trialinfo=cat(1,data_all.trialinfo,TFR_this_session.trialinfo);
                end
            end
        end
    end
end


        
TFR_stimulus           = data_all;
TFR_stimulus.time      = TFR_stimulus.time(25:101);
TFR_stimulus.powspctrm = TFR_stimulus.powspctrm(:,:,:,25:101); 
        
TFR_baseline           = TFR_stimulus;
TFR_baseline.powspctrm = repmat(mean(data_all.powspctrm(:,:,:,11:21),4),1,1,1,size(TFR_baseline.powspctrm,4));
    
save([Path_variables,'/Figure2/Significant_Sensors_Stimulus.mat'],'TFR_stimulus')

save([Path_variables,'/Figure2/Significant_Sensors_Baseline.mat'],'TFR_baseline')

cfg = [];
cfg.channel          = {'MEG'};
cfg.method           = 'montecarlo';
cfg.frequency        = [39 41];
cfg.statistic        = 'ft_statfun_actvsblT';
cfg.correctm         = 'cluster';
cfg.avgovertime      = 'yes';
cfg.avgoverfreq      = 'yes';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.correcttail      = 'prob';
cfg.alpha            = 0.05;
cfg.numrandomization = 10000;
cfg_neighb.method    = 'distance';
cfg.neighbours = ft_prepare_neighbours(cfg_neighb, data_planar);

        
ntrials                       = size(TFR_stimulus.powspctrm,1);
design                        = zeros(2,2*ntrials);
design(1,1:ntrials)           = 1;
design(1,ntrials+1:2*ntrials) = 2;
design(2,1:ntrials)           = [1:ntrials];
design(2,ntrials+1:2*ntrials) = [1:ntrials];
        
cfg.design   = design;
cfg.ivar     = 1;
cfg.uvar     = 2;

%cluster-based statistics
[stat] = ft_freqstatistics(cfg, TFR_stimulus,TFR_baseline);

cd([Path_variables,'/Figure2'])
save([Path_variables,'/Figure2/Significant_Channels_stat.mat'],'stat')