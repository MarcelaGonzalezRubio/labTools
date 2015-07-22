function results = getResultsIncline(SMatrix,params,groups,maxPerturb,plotFlag,indivFlag)

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

results.TMEarlyLate.avg=[];
results.TMEarlyLate.se=[];

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
    TMEarlyLate=[];
    
    
    
    
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
            tmcatchData=adaptData.getParamInCond(params,'catch');
            tmCatch=[tmCatch; smoothedMax(tmcatchData,transientNumPts,stepAsymData)];
            
            % compute OG after as mean values during strides which cause a
            % maximum deviation from zero in STEP LENGTH ASYMMETRY during
            % 'transientNumPts' consecutive strides within first 10 strides
            %             stepAsymData=adaptData.getParamInCond('stepLengthAsym','OG post');
            %             ogafterData=adaptData.getParamInCond(params,'OG post');
            %             ogafter=[ogafter; smoothedMax(ogafterData(1:10,:),transientNumPts,stepAsymData(1:10))];
            
            % compute TM after-effects same as OG after-effect
            stepAsymData=adaptData.getParamInCond('stepLengthAsym','TM post');
            TMafterEarlyData=adaptData.getParamInCond(params,'TM post');
            TMafterEarly=[TMafterEarly; smoothedMax(TMafterEarlyData(1:10,:),transientNumPts,stepAsymData(1:10))];
            
        else
            
            % calculate TM and OG base in same manner as calculating OG after and TM after
            %             OGbaseData=adaptData.getParamInCond(params,'OG base');
            %             OGbase=[OGbase; nanmean(OGbaseData(1:transientNumPts,:))];
            
            TMbaseData=adaptData.getParamInCond(params,'TM base');
            if isempty(TMbaseData)
                TMbaseData=adaptData.getParamInCond(params,{'slow base','fast base'});
            end
            TMbase=[TMbase; nanmean(TMbaseData(1:transientNumPts,:))];
            
            % compute TM post
            if strcmp((adaptData.metaData.conditionName{6}),('TM post'))
                TMafterEarlyData=adaptData.getParamInCond(params,'TM post');
            elseif strcmp((adaptData.metaData.conditionName{9}),('incline post'))
                TMafterEarlyData=adaptData.getParamInCond(params,'incline post');
            end
            TMafterEarly=[TMafterEarly; nanmean(TMafterEarlyData(1:transientNumPts,:))];
        end
        
        if strcmp((adaptData.metaData.conditionName{6}),('TM post'))
            TMafterLateData=adaptData.getParamInCond(params,'TM post');
        elseif strcmp((adaptData.metaData.conditionName{9}),('incline post'))
            TMafterLateData=adaptData.getParamInCond(params,'incline post');
        end
        TMafterLate=[TMafterLate; nanmean(TMafterLateData((end-5)-steadyNumPts+1:(end-5),:))];
        
        % compute TM steady state before catch (mean of first transinetNumPts of last transinetNumPts+5 strides)
        adapt1Data=adaptData.getParamInCond(params,'adaptation');
        EarlyAdapt=[EarlyAdapt; nanmean(adapt1Data(1:transientNumPts,:))];
        
        % compute TM steady state before OG walking (mean of first steadyNumPts of last steadyNumPts+5 strides)
        adapt2Data=adaptData.getParamInCond(params,'adaptation');
        TMsteadyLate=[TMsteadyLate; nanmean(adapt2Data((end-5)-steadyNumPts+1:(end-5),:))];
        
        % compute average adaptation value before the catch
        %         avgAdaptBC=[avgAdaptBC; nanmean(adapt1Data)];
        
        % compute average adaptation of all adaptation walking (both
        % before and after catch)
        adaptAllData=adaptData.getParamInCond(params,{'adaptation'});
        avgAdaptAll=[avgAdaptAll; nanmean(adaptAllData)];
    end
    
    LateEarlyAdapt=[LateEarlyAdapt;TMsteadyLate-EarlyAdapt];
    
    LateSvAdapt=[LateSvAdapt;bsxfun(@minus,TMsteadyLate,TMsteadyLate(:,2))];
    
    TMEarlyLate=[TMEarlyLate;TMafterLate-TMafterEarly];
    
    
    nSubs=length(subjects);
    
    %     results.OGbase.avg(end+1,:)=nanmean(OGbase,1);
    %     results.OGbase.se(end+1,:)=nanstd(OGbase,1)./sqrt(nSubs);
    
    results.TMbase.avg(end+1,:)=nanmean(TMbase,1);
    results.TMbase.se(end+1,:)=nanstd(TMbase,1)./sqrt(nSubs);
    
    results.AvgAdaptAll.avg(end+1,:)=nanmean(avgAdaptAll,1);
    results.AvgAdaptAll.se(end+1,:)=nanstd(avgAdaptAll,1)./sqrt(nSubs);
    
    results.EarlyAdapt.avg(end+1,:)=nanmean(EarlyAdapt,1);
    results.EarlyAdapt.se(end+1,:)=nanstd(EarlyAdapt,1)./sqrt(nSubs);
    
    results.TMsteadyLate.avg(end+1,:)=nanmean(TMsteadyLate,1);
    results.TMsteadyLate.se(end+1,:)=nanstd(TMsteadyLate,1)./sqrt(nSubs);
    
    results.TMafterEarly.avg(end+1,:)=nanmean(TMafterEarly,1);
    results.TMafterEarly.se(end+1,:)=nanstd(TMafterEarly,1)./sqrt(nSubs);
    
    results.TMafterLate.avg(end+1,:)=nanmean(TMafterLate,1);
    results.TMafterLate.se(end+1,:)=nanstd(TMafterLate,1)./sqrt(nSubs);
    
    results.LateEarlyAdapt.avg(end+1,:)=nanmean(LateEarlyAdapt,1);
    results.LateEarlyAdapt.se(end+1,:)=nanstd(LateEarlyAdapt,1)./sqrt(nSubs);
    
    results.LateSvAdapt.avg(end+1,:)=nanmean(LateSvAdapt,1);
    results.LateSvAdapt.se(end+1,:)=nanstd(LateSvAdapt,1)./sqrt(nSubs);
    
    results.TMEarlyLate.avg(end+1,:)=nanmean(TMEarlyLate,1);
    results.TMEarlyLate.se(end+1,:)=nanstd(TMEarlyLate,1)./sqrt(nSubs);
    
    
    if g==1 %This seems ridiculous, but I don't know of another way to do it without making MATLAB mad. The results.(whatever) structure needs to be in this format to make life easier for using SPSS
        for p=1:length(params)
            %             results.OGbase.indiv.(params{p})=[g*ones(nSubs,1) OGbase(:,p)];
            results.TMbase.indiv.(params{p})=[g*ones(nSubs,1) TMbase(:,p)];
            results.AvgAdaptAll.indiv.(params{p})=[g*ones(nSubs,1) avgAdaptAll(:,p)];
            results.EarlyAdapt.indv.(params{p})=[g*ones(nSubs,1) EarlyAdapt(:,p)];
            results.TMsteadyLate.indiv.(params{p})=[g*ones(nSubs,1) TMsteadyLate(:,p)];
            results.TMafterEarly.indiv.(params{p})=[g*ones(nSubs,1) TMafterEarly(:,p)];
            results.TMafterLate.indiv.(params{p})=[g*ones(nSubs,1) TMafterLate(:,p)];
            results.LateEarlyAdapt.indiv.(params{p})=[g*ones(nSubs,1) LateEarlyAdapt(:,p)];
            results.LateSvAdapt.indiv.(params{p})=[g*ones(nSubs,1) LateSvAdapt(:,p)];
            results.TMEarlyLate.indiv.(params{p})=[g*ones(nSubs,1) TMEarlyLate(:,p)];
            
            
            %             results.Transfer.indiv.(params{p})=[g*ones(nSubs,1) transfer(:,p)];
            %             results.Washout.indiv.(params{p})=[g*ones(nSubs,1) washout(:,p)];
            %             results.Transfer2.indiv.(params{p})=[g*ones(nSubs,1) transfer2(:,p)];
            %             results.Washout2.indiv.(params{p})=[g*ones(nSubs,1) washout2(:,p)];
            
            %             results.OGbase.indiv=[g*ones(nSubs,1) OGbase];
            %             results.TMbase.indiv=[g*ones(nSubs,1) TMbase];
            %             results.AvgAdaptBeforeCatch.indiv=[g*ones(nSubs,1) avgAdaptBC];
            %             results.AvgAdaptAll.indiv=[g*ones(nSubs,1) avgAdaptAll];
            %             results.ErrorsOut.indiv=[g*ones(nSubs,1) errorsOut];
            %             results.TMsteadyLate.indiv=[g*ones(nSubs,1) TMsteadyLateBC];
            %             results.catch.indiv=[g*ones(nSubs,1) tmCatch];
            %             results.TMsteadyLate.indiv=[g*ones(nSubs,1) TMsteadyLate];
            %             results.OGafter.indiv=[g*ones(nSubs,1) ogafter];
            %             results.TMafterEarly.indiv=[g*ones(nSubs,1) TMafterEarly];
            %             results.Transfer.indiv=[g*ones(nSubs,1) transfer];
            %             results.Washout.indiv=[g*ones(nSubs,1) washout];
            %             results.Transfer2.indiv=[g*ones(nSubs,1) transfer2];
            %             results.Washout2.indiv=[g*ones(nSubs,1) washout2];
        end
    else
        for p=1:length(params)
            %             results.OGbase.indiv.(params{p})(end+1:end+nSubs,1:2)=[g*ones(nSubs,1) OGbase(:,p)];
            results.TMbase.indiv.(params{p})(end+1:end+nSubs,1:2)=[g*ones(nSubs,1) TMbase(:,p)];
            results.AvgAdaptAll.indiv.(params{p})(end+1:end+nSubs,1:2)=[g*ones(nSubs,1) avgAdaptAll(:,p)];
            results.EarlyAdapt.indv.(params{p})(end+1:end+nSubs,1:2)=[g*ones(nSubs,1) EarlyAdapt(:,p)];
            results.TMsteadyLate.indiv.(params{p})(end+1:end+nSubs,1:2)=[g*ones(nSubs,1) TMsteadyLate(:,p)];
            results.TMafterEarly.indiv.(params{p})(end+1:end+nSubs,1:2)=[g*ones(nSubs,1) TMafterEarly(:,p)];
            results.TMafterLate.indiv.(params{p})(end+1:end+nSubs,1:2)=[g*ones(nSubs,1) TMafterLate(:,p)];
            results.LateEarlyAdapt.indiv.(params{p})(end+1:end+nSubs,1:2)=[g*ones(nSubs,1) LateEarlyAdapt(:,p)];
            results.LateSvAdapt.indiv.(params{p})(end+1:end+nSubs,1:2)=[g*ones(nSubs,1) LateSvAdapt(:,p)];
            results.TMEarlyLate.indiv.(params{p})(end+1:end+nSubs,1:2)=[g*ones(nSubs,1) TMEarlyLate(:,p)];
            %             results.TMsteadyLate.indiv.(params{p})(end+1:end+nSubs,1:2)=[g*ones(nSubs,1) TMsteadyLate(:,p)];
            %             results.OGafter.indiv.(params{p})(end+1:end+nSubs,1:2)=[g*ones(nSubs,1) ogafter(:,p)];
            
            %             results.Transfer.indiv.(params{p})(end+1:end+nSubs,1:2)=[g*ones(nSubs,1) transfer(:,p)];
            %             results.Washout.indiv.(params{p})(end+1:end+nSubs,1:2)=[g*ones(nSubs,1) washout(:,p)];
            %             results.Transfer2.indiv.(params{p})(end+1:end+nSubs,1:2)=[g*ones(nSubs,1) transfer2(:,p)];
            %             results.Washout2.indiv.(params{p})(end+1:end+nSubs,1:2)=[g*ones(nSubs,1) washout2(:,p)];
            
            %             results.OGbase.indiv(end+1:end+nSubs,:)=[g*ones(nSubs,1) OGbase];
            %             results.TMbase.indiv(end+1:end+nSubs,:)=[g*ones(nSubs,1) TMbase];
            %             results.AvgAdaptBeforeCatch.indiv(end+1:end+nSubs,:)=[g*ones(nSubs,1) avgAdaptBC];
            %             results.AvgAdaptAll.indiv(end+1:end+nSubs,:)=[g*ones(nSubs,1) avgAdaptAll];
            %             results.ErrorsOut.indiv(end+1:end+nSubs,:)=[g*ones(nSubs,1) errorsOut];
            %             results.TMsteadyLate.indiv(end+1:end+nSubs,:)=[g*ones(nSubs,1) TMsteadyLateBC];
            %             results.catch.indiv(end+1:end+nSubs,:)=[g*ones(nSubs,1) tmCatch];
            %             results.TMsteadyLate.indiv(end+1:end+nSubs,:)=[g*ones(nSubs,1) TMsteadyLate];
            %             results.OGafter.indiv(end+1:end+nSubs,:)=[g*ones(nSubs,1) ogafter];
            %             results.TMafterEarly.indiv(end+1:end+nSubs,:)=[g*ones(nSubs,1) TMafterEarly];
            %             results.Transfer.indiv(end+1:end+nSubs,:)=[g*ones(nSubs,1) transfer];
            %             results.Washout.indiv(end+1:end+nSubs,:)=[g*ones(nSubs,1) washout];
            %             results.Transfer2.indiv(end+1:end+nSubs,:)=[g*ones(nSubs,1) transfer2];
            %             results.Washout2.indiv(end+1:end+nSubs,:)=[g*ones(nSubs,1) washout2];
        end
    end
end

%plot stuff
if nargin>4 && plotFlag
    
    % FIRST: plot baseline values against catch and transfer
%     epochs={ 'TMafterEarly','TMafterLate','TMEarlyLate'};
% epochs={'TMEarlyLate'};
%     if nargin>5 %I imagine there has to be a better way to do this...
%         %             barGroups(SMatrix,results,groups,{'netContributionNorm2'},epochs,indivFlag)
%         barGroups(SMatrix,results,groups,params,epochs,indivFlag)
%     else
%         barGroups(SMatrix,results,groups,{'netContributionNorm2'},epochs)
%         %             barGroups(SMatrix,results,groups,params,epochs)
%     end
%     
%     % SECOND: plot average adaptation values?
% %     epochs={'EarlyAdapt','TMsteadyLate','LateEarlyAdapt','LateSvAdapt'};
%       epochs={'TMsteadyLate'};
%     if nargin>5
% %         barGroups(SMatrix,results,groups,params,epochs,indivFlag)
%          barGroups(SMatrix,results,groups,{'netContributionNorm2'},epochs,indivFlag)
%     else
%         barGroups(SMatrix,results,groups,params,epochs)
% %           barGroups(SMatrix,results,groups,{'netContributionNorm2'},epochs)
%     end
%     %
%     %       epochs={'Transfer'};
%     %     if nargin>5
%     %         barGroups(SMatrix,results,groups,params,epochs,indivFlag)
%     %     else
%     %         barGroups(SMatrix,results,groups,params,epochs)
%     %     end

   epochs={'TMsteadyLate'};
        if nargin>5
%             barGroups2(SMatrix,results,groups,{'spatialContributionNorm2','stepTimeContributionNorm2'},epochs,indivFlag)
%                    barGroups2(SMatrix,results,groups,{'netContributionNorm2'},epochs,indivFlag)
                barGroups(SMatrix,results,groups,params,epochs,indivFlag)

        else
            barGroups(SMatrix,results,groups,params,epochs)
        end
%         

          epochs={'LateEarlyAdapt'};
        if nargin>5
            barGroups(SMatrix,results,groups,params,epochs,indivFlag)
%             barGroups2(SMatrix,results,groups,{'spatialContributionNorm2','stepTimeContributionNorm2'},epochs,indivFlag)
        else
            barGroups(SMatrix,results,groups,params,epochs)
        end
        
             epochs={'TMafterEarly'};
        if nargin>5
            barGroups(SMatrix,results,groups,params,epochs,indivFlag)
%             barGroups2(SMatrix,results,groups,{'spatialContributionNorm2','stepTimeContributionNorm2'},epochs,indivFlag)
        else
            barGroups(SMatrix,results,groups,params,epochs)
        end
    
end

