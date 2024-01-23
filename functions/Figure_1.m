% exclude : S6 - session 2, S4 - session 2 (BUTTON BROKEN)
clearvars -except Behav_path Path_figures Path_variables

SubjectID={'1' '2' '4' '5' '6' '7' '8' '10' '11' '12' '15' '17' '18' '19' '20' '22' '24','21','25','23'};
Drug=[3 1 1 2 3 2 3; 3 1 3 2 1 2 0; 1 0 2 1 3 2 0; 2 1 1 3 2 3 0; 2 0 3 1 2 3 0; 1 3 2 3 1 2 0; 3 1 1 2 2 3 0; 1 2 2 3 1 3 1; 1 2 2 3 3 1 1; 3 1 2 1 3 2 0; 1 3 3 2 2 1 0; 2 1 3 2 3 1 2; 1 2 3 2 1 3 1; 2 3 1 2 1 3 0; 1 2 3 1 3 2 0; 3 1 2 2 1 3 0; 3 1 1 2 3 2 0; 2 1 1 2 3 3 0; 3 2 1 3 1 2 0;2 3 2 1 1 3 0];
No_button=nan(20,3,3);
Wrong_button=nan(20,3,3);
Hit=nan(20,3,3);
FA=nan(20,3,3);
RT=nan(20,3,3);


for i=1:20
    for d=1:3
        sessions=find(Drug(i,:)==d);
        for s=1:size(sessions,2)
            if s==3
                path=[Behav_path,'/assr_data_extra/']
            else
                path=[Behav_path,'/assr_data/']
            end
            
            listing= dir([path,cell2mat(SubjectID(i)),'/',num2str(sessions(s)),'/'])
            for ii=1:size(listing,1)
                strfind(listing(ii).name,'assr_task');
                if ans==1;
                    load([path,cell2mat(SubjectID(i)),'/',num2str(sessions(s)),'/',listing(ii).name]);
                end
            end
            
             
            Hit(i,d,s)=nansum(session_struct.results.response==2&session_struct.results.tone==2);
            FA(i,d,s)=nansum(session_struct.results.response==2&session_struct.results.tone==1);
            RT(i,d,s)=nanmean(session_struct.results.choice_rt(session_struct.results.response==2&session_struct.results.tone==2));
            
            % in some sessions subjects pressed button 1 instead of button
            % 2
            
            if (nansum(session_struct.results.response==1&session_struct.results.tone==2&isnan(session_struct.results.choice_rt)==0))>=5;
                
                Wrong_button(i,d,s)=1;
                Hit(i,d,s)=nansum(isnan(session_struct.results.choice_rt)==0&session_struct.results.tone==2);
                FA(i,d,s)=nansum(isnan(session_struct.results.choice_rt)==0&session_struct.results.tone==1);
                RT(i,d,s)=nanmean(session_struct.results.choice_rt(isnan(session_struct.results.choice_rt)==0&session_struct.results.tone==2));
            
            end
           
        end
        
        
        
    end
end


FA_all=nanmean(FA,3);
Hit_all=nanmean(Hit./10.*100,3);
RT_all=nanmean(RT,3);


colors =[
    0    162   255;
    255/2 255/2 255/2
    238 34 12]./255;


% plot

figure
subplot(1,3,1)


for Type=1:3
    bar_p=bar(Type,mean(FA_all(:,Type),1),'facecolor',colors(Type,:))
    bar_p.FaceAlpha=0.6;
    hold on
    errorbar(Type,mean(FA_all(:,Type),1),std(FA_all(:,Type))./sqrt(20),'linewidth',3,'color',colors(Type,:))
    
    
    ylim([0 0.2])
    ylabel('False alarm rate (% error)')
    set(gca,'fontsize',15)
    xticks([])
    
    
end

subplot(1,3,2)


for Type=1:3
    bar_p=bar(Type,mean(Hit_all(:,Type),1),'facecolor',colors(Type,:))
    bar_p.FaceAlpha=0.6;
    hold on
    errorbar(Type,mean(Hit_all(:,Type),1),std(Hit_all(:,Type))./sqrt(20),'linewidth',3,'color',colors(Type,:))
    
    hold on
   
    ylim([80 100])
    
    ylabel('Hit rate (% correct)')
    set(gca,'fontsize',15)
    xticks([])
    
    
end

subplot(1,3,3)


for Type=1:3
    bar_p=bar(Type,mean(RT_all(:,Type),1),'facecolor',colors(Type,:))
    bar_p.FaceAlpha=0.6;
    hold on
    errorbar(Type,mean(RT_all(:,Type),1),std(RT_all(:,Type))./sqrt(20),'linewidth',3,'color',colors(Type,:))
  
    ylabel('Reaction time (s)')
    set(gca,'fontsize',15)
    xticks([])
    
    
end



% permutation tests



permutationTest(FA_all(:,1),FA_all(:,2),1000)
permutationTest(FA_all(:,2),FA_all(:,3),1000)
permutationTest(FA_all(:,1),FA_all(:,3),1000)

permutationTest(RT_all(:,1),RT_all(:,2),1000)
permutationTest(RT_all(:,2),RT_all(:,3),1000)
permutationTest(RT_all(:,1),RT_all(:,3),1000)

permutationTest(Hit_all(:,1),Hit_all(:,2),1000)
permutationTest(Hit_all(:,2),Hit_all(:,3),1000)
permutationTest(Hit_all(:,1),Hit_all(:,3),1000)

cd(Path_figures)
saveas(gcf,'Figure1.png')
