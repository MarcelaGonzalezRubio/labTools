function results = getResultsIncline2(SMatrix,params,groups,maxPerturb,plotFlag,indivFlag)

% define number of points to use for calculating values
catchNumPts = 3; %catch
steadyNumPts = 40; %end of adaptation
transientNumPts = 5; %OG and Washout

if nargin<3 || isempty(groups)
    groups=fields(SMatrix);  %default
end
ngroups=length(groups);


% Initialize values to calculate

results.TMbase.avg=[];
results.TMbase.se=[];

results.AvgAdaptAll.avg=[];
results.AvgAdaptAll.se=[];

results.EarlyAdapt.avg=[];
results.EarlyAdapt.se=[];

results.tmCatch.avg=[];
results.tmCatch.se=[];

results.TMsteadyLate.avg=[];
results.TMsteadyLate.se=[];


results.TMafterEarly.avg=[];
results.TMafterEarly.se=[];

results.TMafterLate.avg=[];
results.TMafterLate.se=[];

results.LateEarlyAdapt.avg=[];
results.LateEarlyAdapt.se=[];

results.LateSvAdapt.avg=[];
results.LateSvAdapt.se=[];

results.TMLateEarly.avg=[];
results.TMLateEarly.se=[];

for g=1:ngroups
    
    % get subjects in group
    subjects=SMatrix.(groups{g}).IDs(:,1);
    
    TMbase=[];
    avgAdaptAll=[];
    TMsteadyLate=[];
    EarlyAdapt=[];
    TMafterEarly=[];
    TMafterLate=[];
    LateEarlyAdapt=[];
    LateSvAdapt=[];
    TMLateEarly=[];
    tmCatch=[];
    
    
    
    for s=1:length(subjects)
        % load subject
        load([subjects{s} 'params.mat'])
       
        % remove baseline bias
        adaptData=adaptData.removeBadStrides;
        adaptData=adaptData.removeBias;
        
        
        if nargin>3 && maxPerturb==1
            
            % compute TM and OG base in same manner as calculating OG after and TM after
            
            stepAsymData=adaptData.getParamInCond('stepLengthAsym','TM base');
            TMbaseData=adaptData.getParamInCond(params,'TM base');
            if isempty(TMbaseData)
                stepAsymData=adaptData.getParamInCond('stepLengthAsym',{'slow base','fast base'});
                TMbaseData=adaptData.getParamInCond(params,{'slow base','fast base'});
            end
            TMbase=[TMbase; smoothedMax(TMbaseData(1:10,:),transientNumPts,stepAsymData(1:10))];
            
            % compute catch as mean value during strides which caused a
            % maximum deviation from zero during 'catchNumPts' consecutive
            % strides
            stepAsymData=adaptData.getParamInCond('stepLengthAsym','catch');
             if isempty(stepAsymData)
                  stepAsymData=adaptData.getParamInCond(params,{'TM post','Post-adap','Washout'});
             end
             
             tmcatchData=adaptData.getParamInCond(params,'catch');
%              if isempty(tmcatchData)
%                   tmcatchData=adaptData.getParamInCond(params,{'TM post','Post-adap','Washout'});
%              end
             
             tmCatch=[tmCatch; smoothedMax(tmcatchData,transientNumPts,stepAsymData)];
            
            % compute OG after as mean values during strides which cause a
            % maximum deviation from zero in STEP LENGTH ASYMMETRY during
            % 'transientNumPts' consecutive strides within first 10 strides
            %             stepAsymData=adaptData.getParamInCond('stepLengthAsym','OG post');
            %             ogafterData=adaptData.getParamInCond(params,'OG post');
            %             ogafter=[ogafter; smoothedMax(ogafterData(1:10,:),transientNumPts,stepAsymData(1:10))];
            
            % compute TM after-effects same as OG after-effect
           
            
%             stepAsymData=adaptData.getParamInCond('stepLengthAsym','catch');
%             if isempty(stepAsymData)
            stepAsymData=adaptData.getParamInCond('stepLengthAsym',{'TM post','Post-adap','Washout'});
%             end

%             TMafterEarlyData=adaptData.getParamInCond(params,'catch');
%             if isempty(TMafterEarlyData)
            TMafterEarlyData=adaptData.getParamInCond(params,{'TM post','Post-adap','Washout'});
%             end
             if g==1 && strcmp(adaptData.subData.ID,{'P0011'})
                  TMafterEarlyData=TMafterEarlyData(2:end,:);
             end

            TMafterEarly=[TMafterEarly; smoothedMax(TMafterEarlyData(1:9,:),transientNumPts,stepAsymData(1:9))];
            
        else

            % calculate TM and OG base in same manner as calculating OG after and TM after
            %             OGbaseData=adaptData.getParamInCond(params,'OG base');
            %             OGbase=[OGbase; nanmean(OGbaseData(1:transientNumPts,:))];

            TMbaseData=adaptData.getParamInCond(params,'TM base');
            if isempty(TMbaseData)
                TMbaseData=adaptData.getParamInCond(params,{'slow base','fast base'});
            end
            TMbase=[TMbase; nanmean(TMbaseData(1:transientNumPts,:))];
            
             tmcatchData=adaptData.getParamInCond(params,'catch');
             
             tmCatch=[tmCatch;nanmean(tmcatchData(1:transientNumPts,:))];

            % compute After-Effect on treadmill
%                 TMafterEarlyData=adaptData.getParamInCond(params,'catch');
%                 if isempty(TMafterEarlyData)
                  TMafterEarlyData=adaptData.getParamInCond(params,{'TM post','Post-adap','Washout'});
%                 end
                
                 if group==1 && strcmp(adaptData.subData.ID,{'P0011'})
                      TMafterEarlyData=TMafterEarlyData(2:end,:)
                 end
            TMafterEarly=[TMafterEarly; nanmean(TMafterEarlyData(1:transientNumPts,:))];
        end
           %Compute last steps post-adapt
          TMafterLateData=adaptData.getParamInCond(params,'TM post');
        if isempty(TMafterLateData)
          TMafterLateData=adaptData.getParamInCond(params,{'Post-adap','Washout'});
        end
         TMafterLate=[TMafterLate; nanmean(TMafterLateData((end-5)-steadyNumPts+1:(end-5),:))];
        
        % compute Early adaptation
        
        adapt1Data=adaptData.getParamInCond(params,'adap');
        if strcmp(adaptData.subData.ID,{'P0011'})
          trialNums=adaptData.getTrialsInCond('adapt');
          adapt1Data=adaptData.getParamInTrial(params,trialNums(2));
        end
        EarlyAdapt=[EarlyAdapt; nanmean(adapt1Data(1:transientNumPts,:))];
        
        % compute TM steady state (mean of first steadyNumPts of last steadyNumPts+5 strides)
        adapt2Data=adaptData.getParamInCond(params,'re-adap');
        if isempty(adapt2Data)
            adapt2Data=adaptData.getParamInCond(params,'adap');
%             if strcmp(adaptData.subData.ID,{'P0011'})
%               trialNums=adaptData.getTrialsInCond('adapt');
%               adapt1Data=adaptData.getParamInTrial(params,trialNums(end));
%             end
        end
        TMsteadyLate=[TMsteadyLate; nanmean(adapt2Data((end-5)-steadyNumPts+1:(end-5),:))];
%         TMsteadyLate=[TMafterEarly;tmCatch; nanmean(adapt2Data((end-5)-steadyNumPts+1:(end-5),:))];
        %Compute Avg adaptation
        adaptAllData=adaptData.getParamInCond(params,'adap');
        avgAdaptAll=[avgAdaptAll; nanmean(adaptAllData)];
    end
    
    %Compute Late - Early adaptation
    LateEarlyAdapt=[LateEarlyAdapt;TMsteadyLate-EarlyAdapt];
    
    %Compute LAte - Sv adaptation
    LateSvAdapt=[LateSvAdapt;bsxfun(@minus,TMsteadyLate,TMsteadyLate(:,1))];
    
    %Compute TM late -TM early   
    TMLateEarly=[TMLateEarly;TMafterLate-TMafterEarly];
    
    
    nSubs=length(subjects);
    
    
    results.TMbase.avg(end+1,:)=nanmean(TMbase,1);
    results.TMbase.se(end+1,:)=nanstd(TMbase,1)./sqrt(nSubs);
    
    results.AvgAdaptAll.avg(end+1,:)=nanmean(avgAdaptAll,1);
    results.AvgAdaptAll.se(end+1,:)=nanstd(avgAdaptAll,1)./sqrt(nSubs);
    
    results.EarlyAdapt.avg(end+1,:)=nanmean(EarlyAdapt,1);
    results.EarlyAdapt.se(end+1,:)=nanstd(EarlyAdapt,1)./sqrt(nSubs);
    
    results.tmCatch.avg(end+1,:)=nanmean(tmCatch,1);
    results.tmCatch.se(end+1,:)=nanstd(tmCatch,1)./sqrt(nSubs);
    
    results.TMsteadyLate.avg(end+1,:)=nanmean(TMsteadyLate,1);
    results.TMsteadyLate.se(end+1,:)=nanstd(TMsteadyLate,1)./sqrt(nSubs);
    
    results.TMafterEarly.avg(end+1,:)=nanmean(TMafterEarly,1);
    results.TMafterEarly.se(end+1,:)=nanstd(TMafterEarly,1)./sqrt(nSubs);
%     
    results.TMafterLate.avg(end+1,:)=nanmean(TMafterLate,1);
    results.TMafterLate.se(end+1,:)=nanstd(TMafterLate,1)./sqrt(nSubs);
    
    results.LateEarlyAdapt.avg(end+1,:)=nanmean(LateEarlyAdapt,1);
    results.LateEarlyAdapt.se(end+1,:)=nanstd(LateEarlyAdapt,1)./sqrt(nSubs);
    
    results.LateSvAdapt.avg(end+1,:)=nanmean(LateSvAdapt,1);
    results.LateSvAdapt.se(end+1,:)=nanstd(LateSvAdapt,1)./sqrt(nSubs);
    
    results.TMLateEarly.avg(end+1,:)=nanmean(TMLateEarly,1);
    results.TMLateEarly.se(end+1,:)=nanstd(TMLateEarly,1)./sqrt(nSubs);
%     
    
    if g==1 %This seems ridiculous, but I don't know of another way to do it without making MATLAB mad. The results.(whatever) structure needs to be in this format to make life easier for using SPSS
        for p=1:length(params)
         
            results.TMbase.indiv.(params{p})=[g*ones(nSubs,1) TMbase(:,p)];
            results.AvgAdaptAll.indiv.(params{p})=[g*ones(nSubs,1) avgAdaptAll(:,p)];
            results.EarlyAdapt.indiv.(params{p})=[g*ones(nSubs,1) EarlyAdapt(:,p)];
            results.tmCatch.indiv.(params{p})=[g*ones(nSubs,1) tmCatch(:,p)];     
            results.TMsteadyLate.indiv.(params{p})=[g*ones(nSubs,1) TMsteadyLate(:,p)];
            results.TMafterEarly.indiv.(params{p})=[g*ones(nSubs,1) TMafterEarly(:,p)];
            results.TMafterLate.indiv.(params{p})=[g*ones(nSubs,1) TMafterLate(:,p)];
            results.LateEarlyAdapt.indiv.(params{p})=[g*ones(nSubs,1) LateEarlyAdapt(:,p)];
            results.LateSvAdapt.indiv.(params{p})=[g*ones(nSubs,1) LateSvAdapt(:,p)];
            results.TMLateEarly.indiv.(params{p})=[g*ones(nSubs,1) TMLateEarly(:,p)];
            
        end
    else
        for p=1:length(params)
     
            results.TMbase.indiv.(params{p})(end+1:end+nSubs,1:2)=[g*ones(nSubs,1) TMbase(:,p)];
            results.AvgAdaptAll.indiv.(params{p})(end+1:end+nSubs,1:2)=[g*ones(nSubs,1) avgAdaptAll(:,p)];
            results.EarlyAdapt.indiv.(params{p})(end+1:end+nSubs,1:2)=[g*ones(nSubs,1) EarlyAdapt(:,p)];
            results.tmCatch.indiv.(params{p})(end+1:end+nSubs,1:2)=[g*ones(nSubs,1) tmCatch(:,p)];
            results.TMsteadyLate.indiv.(params{p})(end+1:end+nSubs,1:2)=[g*ones(nSubs,1) TMsteadyLate(:,p)];
            results.TMafterEarly.indiv.(params{p})(end+1:end+nSubs,1:2)=[g*ones(nSubs,1) TMafterEarly(:,p)];
            results.TMafterLate.indiv.(params{p})(end+1:end+nSubs,1:2)=[g*ones(nSubs,1) TMafterLate(:,p)];
            results.LateEarlyAdapt.indiv.(params{p})(end+1:end+nSubs,1:2)=[g*ones(nSubs,1) LateEarlyAdapt(:,p)];
            results.LateSvAdapt.indiv.(params{p})(end+1:end+nSubs,1:2)=[g*ones(nSubs,1) LateSvAdapt(:,p)];
            results.TMLateEarly.indiv.(params{p})(end+1:end+nSubs,1:2)=[g*ones(nSubs,1) TMLateEarly(:,p)];
           
        end
    end
end

%plot stuff
if nargin>4 && plotFlag
    
    % FIRST: plot baseline values against catch and transfer
    %     epochs={'TM'};
    %     if nargin>5 %I imagine there has to be a better way to do this...
    %         barGroups(SMatrix,results,groups,params,epochs,indivFlag)
    %     else
    %         barGroups(SMatrix,results,groups,params,epochs)
    %     end
    
    % SECOND: plot average adaptation values?
%     epochs={'TMafterEarly','tmCatch'};
%     if nargin>5
%         barGroups2(SMatrix,results,groups,params,epochs,indivFlag)
%     else
%         barGroups2(SMatrix,results,groups,params,epochs)
%     end
%     
    epochs={'EarlyAdapt','TMsteadyLate','LateEarlyAdapt','LateSvAdapt'};
    if nargin>5
        barGroups2(SMatrix,results,groups,params,epochs,indivFlag)
    else
        barGroups2(SMatrix,results,groups,params,epochs)
    end
%     %
          epochs={'TMsteadyLate'};
        if nargin>5
%             barGroups2(SMatrix,results,groups,{'spatialContributionNorm2','stepTimeContributionNorm2'},epochs,indivFlag)
%                    barGroups2(SMatrix,results,groups,{'netContributionNorm2'},epochs,indivFlag)
                barGroups2(SMatrix,results,groups,params,epochs,indivFlag)

        else
            barGroups(SMatrix,results,groups,params,epochs)
        end
% %         
% 
%           epochs={'LateEarlyAdapt'};
%         if nargin>5
%             barGroups2(SMatrix,results,groups,params,epochs,indivFlag)
% %             barGroups2(SMatrix,results,groups,{'spatialContributionNorm2','stepTimeContributionNorm2'},epochs,indivFlag)
%         else
%             barGroups(SMatrix,results,groups,params,epochs)
%         end
%         
             epochs={'TMafterEarly'};
        if nargin>5
            barGroups2(SMatrix,results,groups,params,epochs,indivFlag)
%             barGroups2(SMatrix,results,groups,{'spatialContributionNorm2','stepTimeContributionNorm2'},epochs,indivFlag)
        else
            barGroups(SMatrix,results,groups,params,epochs)
        end
    
end

