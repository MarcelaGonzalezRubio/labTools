function [SmoothedData]=AvgTimeCourse_Whole(adaptDataList,params,conditions, A)
% this version was my first attempt to plot the readaptation stuff,
% AND DO SO WITHOUT CROPPING THE DATA
%adaptDataList must be cell array of 'param.mat' file names
%params is cell array of parameters to plot. List with commas to
%plot on separate graphs or with semicolons to plot on same graph.
%conditions is cell array of conditions to plot
%binwidth is the number of data points to average in time
%indivFlag - set to true to plot individual subject time courses
%indivSubs - must be a cell array of 'param.mat' file names that is
%a subset of those in the adaptDataList. Plots specific subjects
%instead of all subjects.

allValues=cell(4, 4);
allValuesC=cell(4, 4);
allValuesALL=cell(4, 4);
Stride2SS=cell(4,4);
whereIS=cell(4,4);
ToAICAnalysis=[];
Divider=cell(4, 4);
%First: see if adaptDataList is a single subject (char), a cell
%array of subject names (one group of subjects), or a cell array of cell arrays of
%subjects names (several groups of subjects), and put all the
%cases into the same format
if isa(adaptDataList,'cell')
    if isa(adaptDataList{1},'cell')
        auxList=adaptDataList;
    else
        auxList{1}=adaptDataList;
    end
elseif isa(adaptDataList,'char')
    auxList{1}={adaptDataList};
end
Ngroups=length(auxList);

%make sure params is a cell array
if isa(params,'char')
    params={params};
end

%check condition input
if nargin>2
    if isa(conditions,'char')
        conditions={conditions};
    end
else
    load(auxList{1}{1})
    conditions=adaptData.metaData.conditionName; %default
end
for c=1:length(conditions)
    cond{c}=conditions{c}(ismember(conditions{c},['A':'Z' 'a':'z' '0':'9'])); %remove non alphanumeric characters
end

if nargin<4
    binwidth=1;
end

if nargin>5 && isa(indivSubs,'cell')
    if ~isa(adaptDataList{1},'cell')
        indivSubs{1}=indivSubs;
    end
elseif nargin>5 && isa(indivSubs,'char')
    indivSubs{1}={indivSubs};
end

%Load data and determine length of conditions
nConds= length(conditions);
s=1;

for group=1:Ngroups
    for subject=1:length(auxList{group})
        %Load subject
        load(auxList{group}{subject});
        
        %normalize contributions based on combined step lengths
        SLf=adaptData.data.getParameter('stepLengthFast');
        SLs=adaptData.data.getParameter('stepLengthSlow');
        Dist=SLf+SLs;
        contLabels={'spatialContribution','stepTimeContribution','velocityContribution','netContribution'};
        [~,dataCols]=isaParameter(adaptData.data,contLabels);
        for c=1:length(contLabels)
            contData=adaptData.data.getParameter(contLabels(c));
            contData=contData./Dist;
            adaptData.data.Data(:,dataCols(c))=contData;
        end
        
        %EDIT: create contribution error values
        vels=adaptData.data.getParameter('stanceSpeedSlow');
        velf=adaptData.data.getParameter('stanceSpeedFast');
        deltaST=adaptData.data.getParameter('stanceTimeDiff');
        velCont=adaptData.data.getParameter('velocityContribution');
        stepCont=adaptData.data.getParameter('stepTimeContribution');
        spatialCont=adaptData.data.getParameter('spatialContribution');
        Tideal=((vels+velf)./2).*deltaST./Dist;
        Sideal=(-velCont)-Tideal;
        [~,dataCols]=isaParameter(adaptData.data,{'Tgoal','Sgoal'});
        adaptData.data.Data(:,dataCols(1))=Tideal-stepCont;
        adaptData.data.Data(:,dataCols(2))=Sideal-spatialCont;
        
        adaptData = adaptData.removeBias; %CJS
        
        for c=1:nConds
            %                         if false
            if strcmpi(conditions{c},'adaptation')  || strcmpi(conditions{c},'TM post') || strcmpi(conditions{c},'re-adaptation')
                trials=adaptData.getTrialsInCond(conditions{c});
                for t=1:length(trials)
                    dataPts=adaptData.getParamInTrial(params,trials(t));
                    nPoints=size(dataPts,1);
                    if nPoints == 0
                        numPts.(cond{c}).(['trial' num2str(t)])(s)=NaN;
                    else
                        numPts.(cond{c}).(['trial' num2str(t)])(s)=nPoints;
                    end
                    for p=1:length(params)
                        %itialize so there are no inconsistant dimensions or out of bounds errors
                        values(group).(params{p}).(cond{c}).(['trial' num2str(t)])(subject,:)=NaN(1,1000); %this assumes that the max number of data points that could exist in a single condition is 1000
                        values(group).(params{p}).(cond{c}).(['trial' num2str(t)])(subject,1:nPoints)=dataPts(:,p);
                    end
                end
                if min(numPts.(cond{c}).(['trial' num2str(t)]))==0
                    bad=find(numPts.(cond{c}).(['trial' num2str(t)])==0);
                    numPts.(cond{c}).(['trial' num2str(t)])=[];
                end
                %~~~~~~~~~~~~
            end
        end
        s=s+1;
    end
end

for group=1:Ngroups
    Xstart=1;
    lineX=0;
    subjects=auxList{group};
    %for c=1:length(conditions)
        if strcmpi(conditions{c},'adaptation') || strcmpi(conditions{c},'re-adaptation') || strcmpi(conditions{c},'TM post')
            trials=[adaptData.getTrialsInCond(conditions{1}) adaptData.getTrialsInCond(conditions{2})];
            %for t=1:length(trials)
            %                 % Min PTS
            %                 [maxPts,loc]=nanmin(numPts.(cond{c}).(['trial' num2str(t)]));
            %                 while maxPts>1.25*nanmax(numPts.(cond{c}).(['trial' num2str(t)])([1:loc-1 loc+1:end]))
            %                     numPts.(cond{c}).(['trial' num2str(t)])(loc)=nanmean(numPts.(cond{c}).(['trial' num2str(t)])([1:loc-1 loc+1:end])); %do not include min in mean
            %                     [maxPts,loc]=nanmin(numPts.(cond{c}).(['trial' num2str(t)]));
            %                 end
            
            %% Max PTS
            [maxPts(t),loc]=nanmax(numPts.(cond{c}).(['trial' num2str(t)]));
            while maxPts>1.25*nanmax(numPts.(cond{c}).(['trial' num2str(t)])([1:loc-1 loc+1:end]))
                numPts.(cond{c}).(['trial' num2str(t)])(loc)=nanmean(numPts.(cond{c}).(['trial' num2str(t)])([1:loc-1 loc+1:end])); %do not include min in mean
                [maxPts(t),loc]=nanmax(numPts.(cond{c}).(['trial' num2str(t)]));
            end
            
            if maxPts==0
                continue
            end
            
            for p=1:length(params)
                
                %                     %                     How you get all the subjects the same length
                %                                         %allValues=[allValues values(group).(params{p}).(cond{c}).(['trial' num2str(t)])(:,1:maxPts-5)]; %% CJS here is where I am taking the adaptation timecourse to fit
                %
                %                                         if strcmpi(conditions{c},'adaptation')
                %                                         allValues{group, p}= [allValues{group, p} values(group).(params{p}).(cond{c}).(['trial' num2str(t)])(:,1:maxPts-5)]; %% CJS here is where I am taking the adaptation timecourse to fit
                %                                         Divider{group, p}=[ Divider{group, p} maxPts-5];
                %                                         elseif strcmpi(conditions{c},'re-adaptation')
                %                                             allValuesC{group, p}= [allValuesC{group, p} values(group).(params{p}).(cond{c}).(['trial' num2str(t)])(:,1:maxPts-5)]; %% CJS here is where I am taking the adaptation timecourse to fit
                %                                             Divider{group, p}=[ Divider{group, p} maxPts-5];
                %                                         end
                
                % %This is to get the REAL end of adaptation for all
                % %subjects
                % for person=1:size(values(group).(params{p}).(cond{c}).(['trial' num2str(t)])(:,end-maxPts(t):end),1)
                %     ender=find(isnan(values(group).(params{p}).(cond{c}).(['trial' num2str(t)])(person,:))==0, 1, 'last');
                %     allValues{group, p}= [allValues{group, p}(:,:); values(group).(params{p}).(cond{c}).(['trial' num2str(t)])(person,ender-(maxPts(t)-1):ender)]; %% CJS here is where I am taking the adaptation timecourse to fit
                % end
                
                %This is to get the all of adaptation for each subject
                maxPts=sum(maxPts);
                t=1;
                for person=1:size(values(group).(params{p}).(cond{c}).(['trial' num2str(t)])(:,end-maxPts:end),1)
                    ender1=find(isnan(values(group).(params{p}).(cond{1}).(['trial' num2str(1)])(person,:))==0, 1, 'last');
                    ender2=find(isnan(values(group).(params{p}).(cond{1}).(['trial' num2str(2)])(person,:))==0, 1, 'last');
                    ender3=find(isnan(values(group).(params{p}).(cond{1}).(['trial' num2str(3)])(person,:))==0, 1, 'last');
                    ender4=find(isnan(values(group).(params{p}).(cond{1}).(['trial' num2str(4)])(person,:))==0, 1, 'last');
                    ender5=find(isnan(values(group).(params{p}).(cond{2}).(['trial' num2str(1)])(person,:))==0, 1, 'last');
                    ender6=find(isnan(values(group).(params{p}).(cond{2}).(['trial' num2str(2)])(person,:))==0, 1, 'last');
                    
                     
                    allValues{group, p}= [allValues{group, p}(:,:); [values(group).(params{p}).(cond{1}).(['trial' num2str(1)])(person,1:ender1) ...
                        values(group).(params{p}).(cond{1}).(['trial' num2str(2)])(person,1:ender2) ...
                        values(group).(params{p}).(cond{1}).(['trial' num2str(3)])(person,1:ender3) ...
                        values(group).(params{p}).(cond{1}).(['trial' num2str(4)])(person,1:ender4) ...
                        nan(1, 700-ender1-ender2-ender3-ender4)]]; %% CJS here is where I am taking the adaptation timecourse to fit
                    Divider{group, p}=[ Divider{group, p}; ender1 ender2 ender3 ender4 ender5 ender6];
                    
                    allValuesC{group, p}= [allValuesC{group, p}(:,:); [values(group).(params{p}).(cond{2}).(['trial' num2str(1)])(person,1:ender5) ...
                        values(group).(params{p}).(cond{2}).(['trial' num2str(2)])(person,1:ender6) ...
                        nan(1, 700-ender5-ender6)]]; %% CJS here is where I am taking the adaptation timecourse to fit
                    
                    allValuesALL{group, p}= [allValuesALL{group, p}(:,:); [values(group).(params{p}).(cond{1}).(['trial' num2str(1)])(person,1:ender1) ...
                        values(group).(params{p}).(cond{1}).(['trial' num2str(2)])(person,1:ender2) ...
                        values(group).(params{p}).(cond{1}).(['trial' num2str(3)])(person,1:ender3) ...
                        values(group).(params{p}).(cond{1}).(['trial' num2str(4)])(person,1:ender4) ...
                        values(group).(params{p}).(cond{2}).(['trial' num2str(1)])(person,1:ender5) ...
                        values(group).(params{p}).(cond{2}).(['trial' num2str(2)])(person,1:ender6) ...
                        nan(1, 950-ender1-ender2-ender3-ender4-ender5-ender6)]];
                    
                    
                end
            end
            %end
        %end
    end
end

% %%If I want to do a variable smooth
% allValues=mean(allValues, 1);
% for qq=1:size(allValues,1)
% [SmoothedData(qq, :)]=bin_data_Variable(allValues(qq,:)',3, 20);
% end
% display('everything is awesome')


% %%%How to calculate strides to ss
figure
tea=1;
milk=1;
for gr=[1:4]
    % Need to visualize to make sure that this is working correctly
    %subplot(2, 2, gr)
    
    for var=[1 2 4]
        
        for qq=1:size(allValues{gr, var},1)
            
            %Smooth the data:
            %allValuesALL{gr, var}(qq,:)=bin_data_Variable(allValuesALL{gr, var}(qq,:)',3, 20); SmoothType='Whole, BWVar, first not before raw min';%, ';
            allValuesALL{gr, var}(qq,:)=bin_dataV1(allValuesALL{gr, var}(qq,:)',20)'; SmoothType='Whole, BW=20, first not before raw min';
            
            %allValuesALL{gr, var}(qq,:)=bin_dataV1(allValuesALL{gr, var}(qq,:)',5)'; SmoothType='Whole, BW=5, consecutive 15 not before min';
            %allValuesALL{gr, var}(qq,:)=[allValuesALL{gr, var}(qq,:)]; SmoothType='Whole, No Smooth, consecutive 5, not before min';
            
            %Here I am using the final steady state that subjects reached
            if gr==1
                ss=A.TMsteady2.indiv.OA(qq, var);
            elseif gr==2
                ss=A.TMsteady2.indiv.OASV(qq, var);
            elseif gr==3
                ss=A.TMsteady2.indiv.YA(qq, var);
            elseif gr==4
                ss=A.TMsteady2.indiv.YASV(qq, var);
            end
            
            % % % %Here I am shifting the SLasym up
            whereIS{gr, var}(qq, :)=find(allValues{gr, var}(qq, :)==nanmin(allValues{gr, var}(qq, 1:50)),1,  'first');%use the non-smoothed data to shift the curves
            if var==4 %SLasym
                %shifter=A.TMsteady2.indiv.OA(qq, 3);
                shifter=.5;
                allValuesALL{gr, var}(qq, :)=allValuesALL{gr, var}(qq, :)+abs(shifter);
                ss=ss+shifter;
            end
            
           t=find(allValuesALL{gr, var}(qq, :)>=ss*.632);
        
% % %             %Here is where I define if there must be a certain number of
% % %             %consecutive stides
% % %             N = 4; % Required number of consecutive numbers following a first one
% % %             x = diff(t)==1;
% % %             f = find([false,x]~=[x,false]);
% % %             g = find(f(2:2:end)-f(1:2:end-1)>=N,1,'first');
% % %             first_t = t(f(2*g-1)); % First t followed by >=N consecutive numbers
% % %             
% % %             if isempty(first_t)
% % %                  first_t=NaN;
% % %             end
% % %             
% % %             nut=2;
% % %             if var==1
% % %             while first_t<whereIS{gr, var}(qq, :)
% % %                 display('ugh')
% % %                 %x = diff(t(whereIS{gr, var}(qq, :):end))==1;
% % %                 %f = find([false,x]~=[x,false]);
% % %                 gALL = find(f(2:2:end)-f(1:2:end-1)>=N);
% % %                 g=gALL(nut);
% % %                 first_t = t(f(2*g-1)); % First t followed by >=N consecutive numbers
% % %                 nut=nut+1;
% % %             end
% % %             end
            
            first_t=t(1);
            knot=2;
            while first_t<=whereIS{gr, var}(qq, :)%5
                first_t=t(knot);
                knot=knot+1;
            end
            
            Stride2SS{gr, var}=[Stride2SS{gr, var} first_t];
            
            if var==4  
 
                if gr == 1
                    %ToAICAnalysis=[ToAICAnalysis; allValuesALL{gr, var}(qq, :)];
                    figure(1)
                    subplot(3, 4, qq)

                    plot([allValuesALL{gr, var}(qq, :)], 'b.-', 'MarkerSize', 25);hold on
                    %plot(t(1), allValuesALL{gr, var}(qq, (t(1))), '.c', 'MarkerSize', 25); hold on
                    plot(whereIS{gr, var}(qq, :):whereIS{gr, var}(qq, :)+9, allValuesALL{gr, var}(qq, whereIS{gr, var}(qq, :):whereIS{gr, var}(qq, :)+9), 'c.', 'MarkerSize', 25);hold on
                    plot(first_t, allValuesALL{gr, var}(qq, first_t), '.r', 'MarkerSize', 25); hold on
                    line([0 900], [ss ss],'Color', 'k', 'LineWidth', 1)
                    line([0 900], [ss*.632 ss*.632],'Color', 'k', 'LineWidth', 1, 'LineStyle',':')
                    
                    line([Divider{gr, var}(qq) Divider{gr, var}(qq)], [-.1 .1], 'Color', 'k', 'LineStyle','--');
                    line([sum(Divider{gr, var}(qq, 1:2)) sum(Divider{gr, var}(qq, 1:2))], [-.1 .1], 'Color', 'k', 'LineStyle','--');
                    line([sum(Divider{gr, var}(qq, 1:3)) sum(Divider{gr, var}(qq, 1:3))], [-.1 .1], 'Color', 'k', 'LineStyle','--');
                    line([sum(Divider{gr, var}(qq, 1:4)) sum(Divider{gr, var}(qq, 1:4))], [-.1 .1], 'Color', 'c', 'LineStyle','--');
                    line([sum(Divider{gr, var}(qq, 1:5)) sum(Divider{gr, var}(qq, 1:5))], [-.1 .1], 'Color', 'c', 'LineStyle','--');
                    line([sum(Divider{gr, var}(qq, 1:6)) sum(Divider{gr, var}(qq, 1:6))], [-.1 .1], 'Color', 'c', 'LineStyle','--');
                    
                    whoIS=cell2mat(adaptDataList{1, gr}(1, qq));
                    title([whoIS(1:end-10) '  Stride to SS = ' num2str(first_t)]);
                    %ylabel('Spatial')
                    ylabel('temporal')
                    xlabel(['OA ' SmoothType])
                    axis tight

                end
                
                
                if gr == 2
                    %ToAICAnalysis=[ToAICAnalysis; allValuesALL{gr, var}(qq, :)];
                    figure(2)
                    subplot(3, 4, qq)
                    plot([allValuesALL{gr, var}(qq, :)], 'b.-', 'MarkerSize', 25);hold on
                    %plot(t(1), allValuesALL{gr, var}(qq, (t(1))), '.c', 'MarkerSize', 25); hold on
                    plot(whereIS{gr, var}(qq, :):whereIS{gr, var}(qq, :)+9, allValuesALL{gr, var}(qq, whereIS{gr, var}(qq, :):whereIS{gr, var}(qq, :)+9), 'c.', 'MarkerSize', 25);hold on
                    plot(first_t, allValuesALL{gr, var}(qq, first_t), '.r', 'MarkerSize', 25); hold on

                    line([0 900], [ss ss],'Color', 'k', 'LineWidth', 1)
                    line([0 900], [ss*.632 ss*.632],'Color', 'k', 'LineWidth', 1, 'LineStyle',':')
                    
                    line([Divider{gr, var}(qq) Divider{gr, var}(qq)], [-.1 .1], 'Color', 'k', 'LineStyle','--');
                    line([sum(Divider{gr, var}(qq, 1:2)) sum(Divider{gr, var}(qq, 1:2))], [-.1 .1], 'Color', 'k', 'LineStyle','--');
                    line([sum(Divider{gr, var}(qq, 1:3)) sum(Divider{gr, var}(qq, 1:3))], [-.1 .1], 'Color', 'k', 'LineStyle','--');
                    line([sum(Divider{gr, var}(qq, 1:4)) sum(Divider{gr, var}(qq, 1:4))], [-.1 .1], 'Color', 'c', 'LineStyle','--');
                    line([sum(Divider{gr, var}(qq, 1:5)) sum(Divider{gr, var}(qq, 1:5))], [-.1 .1], 'Color', 'c', 'LineStyle','--');
                    line([sum(Divider{gr, var}(qq, 1:6)) sum(Divider{gr, var}(qq, 1:6))], [-.1 .1], 'Color', 'c', 'LineStyle','--');
                    
                    whoIS=cell2mat(adaptDataList{1, gr}(1, qq));
                    title([whoIS(1:end-10) '  Stride to SS = ' num2str(first_t)]);
                    %ylabel('Spatial')
                    ylabel('temporal')
                    xlabel(['OASV ' SmoothType])
                end
                
                if gr == 3
                    %ToAICAnalysis=[ToAICAnalysis; allValuesALL{gr, var}(qq, :)];
                    figure(3)
                    subplot(3, 4, qq)
                                        plot([allValuesALL{gr, var}(qq, :)], 'b.-', 'MarkerSize', 25);hold on
                    %plot(t(1), allValuesALL{gr, var}(qq, (t(1))), '.c', 'MarkerSize', 25); hold on
                    plot(whereIS{gr, var}(qq, :):whereIS{gr, var}(qq, :)+9, allValuesALL{gr, var}(qq, whereIS{gr, var}(qq, :):whereIS{gr, var}(qq, :)+9), 'c.', 'MarkerSize', 25);hold on
                    plot(first_t, allValuesALL{gr, var}(qq, first_t), '.r', 'MarkerSize', 25); hold on
                    
                    line([0 900], [ss ss],'Color', 'k', 'LineWidth', 1)
                    line([0 900], [ss*.632 ss*.632],'Color', 'k', 'LineWidth', 1, 'LineStyle',':')
                    
                    line([Divider{gr, var}(qq) Divider{gr, var}(qq)], [-.1 .1], 'Color', 'k', 'LineStyle','--');
                    line([sum(Divider{gr, var}(qq, 1:2)) sum(Divider{gr, var}(qq, 1:2))], [-.1 .1], 'Color', 'k', 'LineStyle','--');
                    line([sum(Divider{gr, var}(qq, 1:3)) sum(Divider{gr, var}(qq, 1:3))], [-.1 .1], 'Color', 'k', 'LineStyle','--');
                    line([sum(Divider{gr, var}(qq, 1:4)) sum(Divider{gr, var}(qq, 1:4))], [-.1 .1], 'Color', 'c', 'LineStyle','--');
                    line([sum(Divider{gr, var}(qq, 1:5)) sum(Divider{gr, var}(qq, 1:5))], [-.1 .1], 'Color', 'c', 'LineStyle','--');
                    line([sum(Divider{gr, var}(qq, 1:6)) sum(Divider{gr, var}(qq, 1:6))], [-.1 .1], 'Color', 'c', 'LineStyle','--');
                    
                    whoIS=cell2mat(adaptDataList{1, gr}(1, qq));
                    title([whoIS(1:end-10) '  Stride to SS = ' num2str(first_t)]);
                    ylabel('Spatial')
                    xlabel(['YA ' SmoothType])
                end
                
                
                if gr == 4
                    ToAICAnalysis=[ToAICAnalysis; allValuesALL{gr, var}(qq, :)];
                    figure(4)
                    subplot(3, 4, qq)
                    
                    plot([allValuesALL{gr, var}(qq, :)], 'b.-', 'MarkerSize', 25);hold on
                    %plot(t(1), allValuesALL{gr, var}(qq, (t(1))), '.c', 'MarkerSize', 25); hold on
                    plot(whereIS{gr, var}(qq, :):whereIS{gr, var}(qq, :)+9, allValuesALL{gr, var}(qq, whereIS{gr, var}(qq, :):whereIS{gr, var}(qq, :)+9), 'c.', 'MarkerSize', 25);hold on
                    plot(first_t, allValuesALL{gr, var}(qq, first_t), '.r', 'MarkerSize', 25); hold on

                    line([0 900], [ss ss],'Color', 'k', 'LineWidth', 1)
                    line([0 900], [ss*.632 ss*.632],'Color', 'k', 'LineWidth', 1, 'LineStyle',':')
                    
                    line([Divider{gr, var}(qq) Divider{gr, var}(qq)], [-.1 .1], 'Color', 'k', 'LineStyle','--');
                    line([sum(Divider{gr, var}(qq, 1:2)) sum(Divider{gr, var}(qq, 1:2))], [-.1 .1], 'Color', 'k', 'LineStyle','--');
                    line([sum(Divider{gr, var}(qq, 1:3)) sum(Divider{gr, var}(qq, 1:3))], [-.1 .1], 'Color', 'k', 'LineStyle','--');
                    line([sum(Divider{gr, var}(qq, 1:4)) sum(Divider{gr, var}(qq, 1:4))], [-.1 .1], 'Color', 'c', 'LineStyle','--');
                    line([sum(Divider{gr, var}(qq, 1:5)) sum(Divider{gr, var}(qq, 1:5))], [-.1 .1], 'Color', 'c', 'LineStyle','--');
                    line([sum(Divider{gr, var}(qq, 1:6)) sum(Divider{gr, var}(qq, 1:6))], [-.1 .1], 'Color', 'c', 'LineStyle','--');
                    
                    whoIS=cell2mat(adaptDataList{1, gr}(1, qq));
                    title([whoIS(1:end-10) '  Stride to SS = ' num2str(first_t)]);
                    ylabel('Spatial')
                    xlabel(['YASV ' SmoothType])
                end
            end          
        end
        tea=1;
        milk=1;
    end
    %     subplot(3, 4, qq+1)
    %     ylabel('Spatial')
    %     title('YASV, Variable Smoothed, First Stride ')
    %legend('Data', 'First Stride', 'First Stride (Ignore first 5)')
    %title(num2str(gr))
end

               figure(5)
               plot(mean(allValuesALL{1, 1}), '.', 'MarkerSize', 25);hold on
                plot(mean(allValuesALL{2, 1}), '.', 'MarkerSize', 25);hold on
                plot(mean(allValuesALL{3, 1}), '.', 'MarkerSize', 25);hold on
                plot(mean(allValuesALL{4, 1}), '.', 'MarkerSize', 25);hold on
                legend('OA', 'OASV', 'YA', 'YASV')
                 line([Divider{gr, var}(1) Divider{gr, var}(1)], [-.1 .25], 'Color', 'k', 'LineStyle','--');
                    line([sum(Divider{gr, var}(1:2)) sum(Divider{gr, var}(1:2))], [-.1 .25], 'Color', 'k', 'LineStyle','--');
                    line([sum(Divider{gr, var}(1:3)) sum(Divider{gr, var}(1:3))], [-.1 .25], 'Color', 'k', 'LineStyle','--');
                    line([sum(Divider{gr, var}(1:4)) sum(Divider{gr, var}(1:4))], [-.1 .25], 'Color', 'c', 'LineStyle','--');
                    line([sum(Divider{gr, var}(1:5)) sum(Divider{gr, var}(1:5))], [-.1 .25], 'Color', 'c', 'LineStyle','--');
                    line([sum(Divider{gr, var}(1:6)) sum(Divider{gr, var}(1:6))], [-.1 .25], 'Color', 'c', 'LineStyle','--');
                    title(['Group Averages: ' SmoothType]);
                    ylabel('Spatial')
                    xlabel(['Strides'])

ReadyAIC=mean(ToAICAnalysis);
display('Milk and Honey Dude')

% % %%%How to start the subjects at the same place
% for gr=[1:4]
%     for var=[1:4]
%         for qq=1:size(allValues{gr, var},1)
%             allValues{gr, var}(qq, :)= [allValues{gr, var}(qq, :)-mean(allValues{gr, var}(qq, 1:3))]; %% CJS here is where I am taking the adaptation timecourse to fit
%         end
%     end
% end

% figure%Just to look at OASV
% for i=1:size(allValues{2, 1},1)
%     plot(allValues{2, 1}(i, 1:30)); hold on
% end
% plot(mean(allValues{2, 1}(:,1:30)), 'k', 'Linewidth', 5)
% title('OASV Individual Subjects')
% ylabel('Spatial')

% figure
% plot(mean(allValues{1, 1}(:,:)), '.', 'MarkerSize', 20); hold on
% plot(mean(allValues{2, 1}(:,:)), '.', 'MarkerSize', 20); hold on
% plot(mean(allValues{3, 1}(:,:)), '.', 'MarkerSize', 20); hold on
% plot(mean(allValues{4, 1}(:,:)), '.', 'MarkerSize', 20); hold on
% % plot(mean(allValues{1, 1}(:,1:30)), '.', 'MarkerSize', 20); hold on
% % plot(mean(allValues{2, 1}(:,1:30)), '.', 'MarkerSize', 20); hold on
% % plot(mean(allValues{3, 1}(:,1:30)), '.', 'MarkerSize', 20); hold on
% % plot(mean(allValues{4, 1}(:,1:30)), '.', 'MarkerSize', 20); hold on
% ylabel('Spatial (Not Shifted)')
% legend('OA', 'OASV', 'YA', 'YASV')

% %%If I want to do a variable smooth
% %allValues=mean(allValues, 1);
% [SmoothedData(1, :)]=bin_dataV1(mean(allValues{1, 1})',3);
% [SmoothedData(2, :)]=bin_dataV1(mean(allValues{2, 1})',3);
% [SmoothedData(3, :)]=bin_dataV1(mean(allValues{3, 1})',3);
% [SmoothedData(4, :)]=bin_dataV1(mean(allValues{4, 1})',3);
% display('everything is awesome')
%
% figure
% plot(SmoothedData(1, :)); hold on
% plot(SmoothedData(2, :)); hold on
% plot(SmoothedData(3, :)); hold on
% plot(SmoothedData(4, :)); hold on
% ylabel('Spatial (Smoothed)')
%  legend('OA', 'OASV', 'YA', 'YASV')

% % % %Get ready for stats by binning data
% % % NumBin=1;
% % % BWnum=30;
% % % SoYeah=cell(4,4);
% % % %for ee=1:floor(size(allValues,2)/30)
% % % for pa=1:4
% % %     for g=1:4
% % %         %for ee=[0:30:30*floor(size(allValues{g, pa},2)/30)-1]
% % %         for ee=1:NumBin%floor(size(allValues{g, pa},2)/30)
% % %             for pp=1:size(allValues{g, pa},1)
% % %                 %%%%%SoYeah{g, pa}=[SoYeah; nanmean(allValues{g, pa}(pp,ee+1:ee+29))];
% % %                 %%%%%SoYeah(pp, ee)=nanmean(allValues(pp,ee*30-29:ee*30));
% % %                 SoYeah{g, pa}(pp, ee)=nanmean(allValues{g, pa}(pp,(BWnum*ee)-(BWnum-1):BWnum*ee));%REAL
% % %                 %                 %If I want to access the SS1 when things are shifted
% % %                 %                 SoYeah{g, pa}(pp, ee)=nanmean(allValues{g, pa}(pp,(end-5)-(40)+1:(end-5),:));%If I want to access the SS1 when things are shifted
% % %             end
% % %         end
% % %     end
% % % end
% % % 
IntoStata_dataS= [];
IntoStata_dataT= [];
IntoStata_dataV= [];
IntoStata_dataN= [];
IntoStata_epoch= [];
IntoStata__group=[];
IntoStata_p=[];
IntoStata_visit=[];
IntoStata_age=[];
IntoStata_strideS=[];
IntoStata_strideT=[];
IntoStata_strideV=[];
IntoStata_strideN=[];
IntoStata_group=[];
% % % 
for g=1:size(adaptDataList, 2)%4
% % %     for ee=1:NumBin%floor(size(allValues{g, pa},2)/30)
% % %         IntoStata_dataS= [IntoStata_dataS; SoYeah{g, 1}(:,ee)];
% % %         IntoStata_dataT= [IntoStata_dataT; SoYeah{g, 2}(:,ee)];
% % %         IntoStata_dataV= [IntoStata_dataV; SoYeah{g, 3}(:,ee)];
% % %         IntoStata_dataN= [IntoStata_dataN; SoYeah{g, 4}(:,ee)];
% % %         IntoStata_epoch= [IntoStata_epoch; ee*ones(size(SoYeah{g, 1}(:,ee)))];
% % %         
        IntoStata_strideS= [IntoStata_strideS; Stride2SS{g, 1}'];
        IntoStata_strideT= [IntoStata_strideT; Stride2SS{g, 2}'];
        IntoStata_strideV= [IntoStata_strideV; Stride2SS{g, 3}'];
        IntoStata_strideN= [IntoStata_strideN; Stride2SS{g, 4}'];
        
        
        
        
        %             if g==1  %ISU
        %                 IntoStata_group=[IntoStata_group;  1*ones(size(SoYeah{g, 1}(:,ee)))];
        %
        %             elseif g==2 %ISD
        % IntoStata_group=[IntoStata_group;  2*ones(size(SoYeah{g, 1}(:,ee)))];
        %
        %             elseif g==3 %ISF
        % IntoStata_group=[IntoStata_group;  3*ones(size(SoYeah{g, 1}(:,ee)))];
        %
        %             end
% % %         
% % %         if g<=2 %Old person
% % %             IntoStata_age=[IntoStata_age;  1*ones(size(SoYeah{g, 1}(:,ee)))];
% % %             if g==1  %OA
% % %                 IntoStata_p=[IntoStata_p; [1:11]' ];
% % %                 IntoStata_visit=[IntoStata_visit;  1*ones(size(SoYeah{g, 1}(:,ee)))];
% % %                 IntoStata_group=[ IntoStata_group; 1*ones(size(SoYeah{g, 1}(:,ee)))];
% % %             elseif g==2 %OASV
% % %                 IntoStata_p=[IntoStata_p; [1:8]' ];
% % %                 IntoStata_visit=[IntoStata_visit;  2*ones(size(SoYeah{g, 1}(:,ee)))];
% % %                 IntoStata_group=[ IntoStata_group; 2*ones(size(SoYeah{g, 1}(:,ee)))];
% % %             end
% % %             
% % %         else %Young Person
% % %             IntoStata_age=[IntoStata_age; 0*ones(size(SoYeah{g, 1}(:,ee)))];
% % %             if g==3 %YA
% % %                 IntoStata_p=[IntoStata_p; [12:19]' ];
% % %                 IntoStata_visit=[IntoStata_visit;  1*ones(size(SoYeah{g, 1}(:,ee)))];
% % %                 IntoStata_group=[ IntoStata_group; 3*ones(size(SoYeah{g, 1}(:,ee)))];
% % %             elseif g==4 %YASV
% % %                 IntoStata_p=[IntoStata_p; 13; 15; 16; 17; 18; 19];
% % %                 IntoStata_visit=[IntoStata_visit;  2*ones(size(SoYeah{g, 1}(:,ee)))];
% % %                 IntoStata_group=[ IntoStata_group; 4*ones(size(SoYeah{g, 1}(:,ee)))];
% % %             end
% % %         end
    end
%%%end

FinalStata= [IntoStata_epoch IntoStata_p IntoStata_age IntoStata_visit IntoStata_dataS IntoStata_dataT IntoStata_dataV IntoStata_dataN IntoStata_group IntoStata_strideS IntoStata_strideT ones(33,1) IntoStata_strideN];

display('everything is awesome')

% key=1;
%
% figure
% plot(allValues(key,:), '.k', 'MarkerSize', 25)
% hold on
% plot(SmoothedData(key,:), '.r', 'MarkerSize', 25)
% legend('Raw', 'Variable Smooth')
end

