function results = barGroupsFirstSteps(SMatrix,params,groups,indivFlag)

% Set colors
poster_colors;
% Set colors order
ColorOrder=[p_plum; p_orange; p_fade_green; p_fade_blue; p_red; p_green; p_blue; p_fade_red; p_lime; p_yellow; p_gray; p_black; p_red];
%ColorOrder=[p_red; p_orange; p_green; p_blue;   p_fade_green; p_fade_blue;  p_fade_red; p_yellow;p_blue; p_orange; p_green;p_red; p_orange; p_green; p_blue;   p_fade_green; p_fade_blue;  p_fade_red; p_yellow;p_blue; p_orange; p_green];
             
catchNumPts = 3; % catch
steadyNumPts = 40; %end of adaptation
transientNumPts = 5; % OG and Washout

if nargin<3 || isempty(groups)
    groups=fields(SMatrix);          
end
ngroups=length(groups);

results.TMstart.avg=[];
results.TMstart.sd=[];
results.TMsteady1.avg=[];
results.TMsteady1.sd=[];

results.PerForget.avg=[];
results.PerForget.sd=[];
results.catch.avg=[];
results.catch.sd=[];
results.TMsteady2.avg=[];
results.TMsteady2.sd=[];
% results.Strides2SS.avg=[171.1818182	35.54545455	NaN	135; 196.125	26.75	NaN	87.125; 69.375	37.875	NaN	92.75; 66	62.33333333	NaN	82];
% results.Strides2SS.sd=[28.9839511	19.75824966	NaN	19.86957473; 34.90161555	21.50975692	NaN	21.07591508; 27.36390839	17.51166703	NaN	19.91656704; 40.44584857	36.0172798	NaN	30.19050623];
% results.Strides2SS.avg=[171.1818182	35.54545455	1	135; 196.125	26.75	1	87.125; 69.375	37.875	1	92.75 ;66	62.33333333	1	82];
% results.Strides2SS.sd=[28.9839511	19.75824966	0	19.86957473; 34.90161555	21.50975692	0	21.07591508; 27.36390839	17.51166703	0	19.91656704; 40.44584857	36.0172798	0	30.19050623];
% results.Strides2SS.avg=[107.7272727	22	1	95.81818182; 45.625	26.125	1	49.625; 39.375	20.625	1	76.125; 34.83333333	37.5	1	42.83333333];
% results.Strides2SS.sd=[21.86865374	9.654956334	0	15.7559906; 12.58391476	14.62194718	0	12.88816941; 17.83349327	7.096119412	0	22.67821288; 19.96900376	28.36400301	0	18.53899794];
results.Strides2SS.avg=[149.1818182	17.63636364	1	126.0909091; 141.125	27	1	64.5; 52.71428571	32.85714286	1	40.85714286; 17.2	26.2	1	15.8];
results.Strides2SS.sd=[28.90102804	5.233237703	0	24.21661884; 28.51530761	14.04584331	0	17.32772345; 20.53882999	8.327889378	0	15.2961691; 5.180733539	11.69786305	0	6.923871749];
results.Strides2SS.indiv.(groups{1})=[52	32	1	48; ...
188	55	1	189; ...
225	2	1	244; ...
186	14	1	185; ...
169	2	1	159; ...
124	9	1	96; ...
12	8	1	18; ...
200	15	1	87; ...
54	9	1	39; ...
345	7	1	238; ...
86	41	1	84];
results.Strides2SS.indiv.(groups{2})=[129	6	1	32; ...
194	122	1	129; ...
261	9	1	114; ...
160	6	1	25; ...
125	31	1	27; ...
56	3	1	60; ...
9	12	1	9; ...
195	27	1	120];
results.Strides2SS.indiv.(groups{3})=[2	46	1	2
51	50	1	29
153	29	1	125
97	18	1	54
33	10	1	33
24	10	1	29
9	67	1	14];
results.Strides2SS.indiv.(groups{4})=[3	4	1	2; ...
17	27	1	3; ...
33	4	1	33; ...; ...
23	68	1	32; ...
10	28	1	9];

results.OGafter.avg=[];
results.OGafter.sd=[];
results.TMafter.avg=[];
results.TMafter.sd=[];
results.Transfer.avg=[];
results.Transfer.sd=[];
results.Washout.avg=[];
results.Washout.sd=[];
results.Transfer2.avg=[];
results.Transfer2.sd=[];
results.Washout2.avg=[];
results.Washout2.sd=[];
results.Remember.avg=[];
results.Remember.sd=[];

results.MagAdapt1.avg=[];
results.MagAdapt1.sd=[];
results.Forget.avg=[];
results.Forget.sd=[];
results.AVGForget.avg=[];
results.AVGForget.sd=[];

for g=1:ngroups
    %get subjects in group
    subjects=SMatrix.(groups{g}).IDs(:,1);
    
    TMstart=[];
    perALL=[];
    remember=[];
    stepsymmetryCatch=[];
    perforget=[];
    AVGforget=[];
    forget=[];
    tmsteady1=[];
    tmcatch=[];
    tmsteady2=[];
    ogafter=[];
    tmafter=[];
    transfer=[];
    washout=[];
    transfer2=[];
    washout2=[];
    MagAdapt1=[];
        
    for s=1:length(subjects)
        %load subject
        load([subjects{s} 'params.mat'])
     
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
                
        %remove baseline bias
        adaptData=adaptData.removeBias;
        
        %Starting Adaptation
        TMstartData=adaptData.getParamInCond(params,'adaptation');
        TMstart=[TMstart;nanmean(TMstartData(1:5,:), 1)];
        
%          %calculate TM steady state at the end of block 3
%        Atrial=getTrialsInCond(adaptData,'adaptation')
%         tmsteadyTRIALSPECIFIC=[adaptData.getParamInTrial(params(1),Atrial(3)) adaptData.getParamInTrial(params(2),Atrial(3)) adaptData.getParamInTrial(params(3),Atrial(3)) adaptData.getParamInTrial(params(4),Atrial(3))];
%         tmsteady1=[tmsteady1;nanmean(tmsteadyTRIALSPECIFIC((end-5)-(steadyNumPts)+1:(end-5),:))];
        
%         %calculate TM steady state #1
         tmsteady1Data=adaptData.getParamInCond(params,'adaptation');
         tmsteady1=[tmsteady1;nanmean(tmsteady1Data((end-5)-(steadyNumPts)+1:(end-5),:))];
        
       
        %calculate forgetting between B1 and B2 of adaptaiton
        test=adaptData.metaData.conditionName;
        test(cellfun(@isempty,test))={''};
        epoch=find(ismember(test, 'adaptation')==1);
        wantedtrials=adaptData.metaData.trialsInCondition{epoch};  
        forgetB1Data=adaptData.getParamInTrial(params,wantedtrials(1));
        forgetB2Data=adaptData.getParamInTrial(params,wantedtrials(2));
        forget=[forget; nanmean(forgetB1Data(end-4:end,:))-nanmean(forgetB2Data(1:5,:))];
        
        %calculate the average forgetting between adaptation and
        %re-adaptation
        forgetB3Data=adaptData.getParamInTrial(params,wantedtrials(3));
        forgetB4Data=adaptData.getParamInTrial(params,wantedtrials(4));
        BOCHI=[nanmean(forgetB1Data(end-4:end,:))-nanmean(forgetB2Data(1:5,:));...
            nanmean(forgetB2Data(end-4:end,:))-nanmean(forgetB3Data(1:5,:));...
            nanmean(forgetB3Data(end-4:end,:))-nanmean(forgetB4Data(1:5,:))];
        AVGforget=[AVGforget; nanmean(BOCHI)];
       
%         per=[(nanmean(forgetB1Data(end-4:end,:))-nanmean(forgetB2Data(1:5,:)))./nanmean(forgetB1Data(end-4:end,:));...
%             (nanmean(forgetB2Data(end-4:end,:))-nanmean(forgetB3Data(1:5,:)))./nanmean(forgetB2Data(end-4:end,:));...
%             (nanmean(forgetB3Data(end-4:end,:))-nanmean(forgetB4Data(1:5,:)))./nanmean(forgetB3Data(end-4:end,:))];

%    per=[(nanmean(forgetB1Data(end-10:end-5,:))-nanmean(forgetB2Data(1:5,:)))./nanmean(forgetB1Data(end-10:end-5,:));...
%             (nanmean(forgetB2Data(end-10:end-5,:))-nanmean(forgetB3Data(1:5,:)))./nanmean(forgetB2Data(end-10:end-5,:));...
%             (nanmean(forgetB3Data(end-10:end-5,:))-nanmean(forgetB4Data(1:5,:)))./nanmean(forgetB3Data(end-10:end-5,:))];

% % ***** 20-10 and 5-3
%            per=[(nanmean(forgetB1Data(end-29:end-10,:))-nanmean(forgetB2Data(4:8,:)))./nanmean(forgetB1Data(end-29:end-10,:));...
%             (nanmean(forgetB2Data(end-29:end-10,:))-nanmean(forgetB3Data(4:8,:)))./nanmean(forgetB2Data(end-29:end-10,:));...
%             (nanmean(forgetB3Data(end-29:end-10,:))-nanmean(forgetB4Data(4:8,:)))./nanmean(forgetB3Data(end-29:end-10,:))];

% % ***** Add Minimum, only for the net 
% minValue=abs(nanmin([forgetB1Data; forgetB2Data; forgetB3Data; forgetB4Data]));
% minValue(1:3)=0;
% forgetB1Data=forgetB1Data+repmat(minValue,length(forgetB1Data),1);
% forgetB2Data=forgetB2Data+repmat(minValue,length(forgetB2Data),1);
% forgetB3Data=forgetB3Data+repmat(minValue,length(forgetB3Data),1);
% forgetB4Data=forgetB4Data+repmat(minValue,length(forgetB4Data),1);

% ***** Add constant, only for the net 04/2015
minValue=[ 0 0 0 0.4];
forgetB1Data=forgetB1Data+repmat(minValue,length(forgetB1Data),1);
forgetB2Data=forgetB2Data+repmat(minValue,length(forgetB2Data),1);
forgetB3Data=forgetB3Data+repmat(minValue,length(forgetB3Data),1);
forgetB4Data=forgetB4Data+repmat(minValue,length(forgetB4Data),1);

           per=[(nanmean(forgetB1Data(end-29:end-10,:))-nanmean(forgetB2Data(4:8,:)))./nanmean(forgetB1Data(end-29:end-10,:));...
            (nanmean(forgetB2Data(end-29:end-10,:))-nanmean(forgetB3Data(4:8,:)))./nanmean(forgetB2Data(end-29:end-10,:));...
            (nanmean(forgetB3Data(end-29:end-10,:))-nanmean(forgetB4Data(4:8,:)))./nanmean(forgetB3Data(end-29:end-10,:))];



% % ***** Beggining Ratios
%            per=[(nanmean(forgetB2Data(4:8,:)))./nanmean(forgetB1Data(4:8,:));...
%             (nanmean(forgetB3Data(4:8,:)))./nanmean(forgetB2Data(4:8,:));...
%             (nanmean(forgetB4Data(4:8,:)))./nanmean(forgetB3Data(4:8,:))];

% ***** End Ratios        
% % % %             per=[(nanmean(forgetB2Data(11:15,:)))./nanmean(forgetB1Data(end-29:end-10,:));...
% % % %             (nanmean(forgetB3Data(11:15,:)))./nanmean(forgetB2Data(end-29:end-10,:));...
% % % %             (nanmean(forgetB4Data(11:15,:)))./nanmean(forgetB3Data(end-29:end-10,:))];


%             per=[(nanmean(forgetB2Data(4:8,:)))./nanmean(forgetB1Data(end-29:end-10,:));...
%             (nanmean(forgetB3Data(4:8,:)))./nanmean(forgetB2Data(end-29:end-10,:));...
%             (nanmean(forgetB4Data(4:8,:)))./nanmean(forgetB3Data(end-29:end-10,:))];

% ***** Long Shot       
%            per=[(nanmean(forgetB1Data(end-14:end-10,:))-nanmean(forgetB2Data(4:8,:)))./(nanmean(forgetB1Data(end-14:end-10,:))-nanmean(forgetB1Data(4:8,:)));...
%             (nanmean(forgetB2Data(end-14:end-10,:))-nanmean(forgetB3Data(4:8,:)))./(nanmean(forgetB2Data(end-14:end-10,:))-nanmean(forgetB2Data(4:8,:)));...
%             (nanmean(forgetB3Data(end-14:end-10,:))-nanmean(forgetB4Data(4:8,:)))./(nanmean(forgetB3Data(end-14:end-10,:))-nanmean(forgetB3Data(4:8,:)))];
% for rat=1:size(params,2)           
% per(:,rat)=[pdist([nanmean(forgetB1Data(end-29:end-10,rat));nanmean(forgetB2Data(4:8,rat))],'euclidean')./pdist([nanmean(forgetB1Data(end-29:end-10,rat));nanmean(forgetB1Data(4:8,rat))],'euclidean');...
%             pdist([nanmean(forgetB2Data(end-29:end-10,rat));nanmean(forgetB3Data(4:8,rat))],'euclidean')./pdist([nanmean(forgetB2Data(end-29:end-10,rat));nanmean(forgetB2Data(4:8,rat))],'euclidean');...
%             pdist([nanmean(forgetB3Data(end-29:end-10,rat));nanmean(forgetB4Data(4:8,rat))],'euclidean')./pdist([nanmean(forgetB3Data(end-29:end-10,rat));nanmean(forgetB3Data(4:8,rat))],'euclidean')];
%           %end-14:end-10
%           %per(:,rat)=[pdist([nanmean(forgetB1Data(end-29:end-10,rat));nanmean(forgetB2Data(4:8,rat))],'euclidean')./pdist([nanmean(forgetB1Data(end-29:end-10,rat));nanmean(forgetB1Data(4:8,rat))],'euclidean')];
% end
        
% per=[sqrt((nanmean(forgetB1Data(end-4:end,:))-nanmean(forgetB2Data(1:5,:))).^2)./nanmean(forgetB1Data(end-4:end,:));...
%             sqrt((nanmean(forgetB2Data(end-4:end,:))-nanmean(forgetB3Data(1:5,:))).^2)./nanmean(forgetB2Data(end-4:end,:));...
%             sqrt((nanmean(forgetB3Data(end-4:end,:))-nanmean(forgetB4Data(1:5,:))).^2)./nanmean(forgetB3Data(end-4:end,:))];
        
        
%         per=[sqrt((nanmean(forgetB1Data(end-4:end,:))-nanmean(forgetB2Data(6:10,:))).^2)./nanmean(forgetB1Data(end-4:end,:));...
%             sqrt((nanmean(forgetB2Data(end-4:end,:))-nanmean(forgetB3Data(6:10,:))).^2)./nanmean(forgetB2Data(end-4:end,:));...
%             sqrt((nanmean(forgetB3Data(end-4:end,:))-nanmean(forgetB4Data(6:10,:))).^2)./nanmean(forgetB3Data(end-4:end,:))];
        
   
%         per=[(nanmean(forgetB1Data(end-4:end,:))-nanmean(forgetB2Data(6:10,:)))./nanmean(forgetB1Data(end-4:end,:));...
%             (nanmean(forgetB2Data(end-4:end,:))-nanmean(forgetB3Data(6:10,:)))./nanmean(forgetB2Data(end-4:end,:));...
%             (nanmean(forgetB3Data(end-4:end,:))-nanmean(forgetB4Data(6:10,:)))./nanmean(forgetB3Data(end-4:end,:))];


%         per=[(nanmean(forgetB1Data(end-9:end-5,:))-nanmean(forgetB2Data(6:10,:)))./nanmean(forgetB1Data(end-9:end-5,:));...
%             (nanmean(forgetB2Data(end-9:end-5,:))-nanmean(forgetB3Data(6:10,:)))./nanmean(forgetB2Data(end-9:end-5,:));...
%             (nanmean(forgetB3Data(end-9:end-5,:))-nanmean(forgetB4Data(6:10,:)))./nanmean(forgetB3Data(end-9:end-5,:))];

perALL=[perALL; per(:,1)];

%per=real(per);
 %perforget=[perforget; 100-(100*(nanmean(per)))];%This is new 4/2015
  perforget=[perforget; (100*(nanmean(per)))];
%perforget=[perforget; 100*((per))];
     


%         %calculate catch as the first 3 strides of adaptation
%         tmcatchData=adaptData.getParamInCond(params,'catch');
%         if isempty(tmcatchData)
%             newtmcatchData=NaN(1,length(params));
%         elseif size(tmcatchData,1)<3
%             newtmcatchData=nanmean(tmcatchData);
%         else
%             newtmcatchData=nanmean(tmcatchData(1:catchNumPts, :));
%         end
%         [~,maxLoc]=max(abs(newtmcatchData),[],1);
%         ind=sub2ind(size(newtmcatchData),maxLoc,1:length(params));
%         tmcatch=[tmcatch; newtmcatchData(ind)];


       %%% ~~~~~~~~
%             calculate catch as mean value during strides which caused a
%             maximum deviation from zero in step length asymmetry during 
%             'catchNumPts' consecutive steps
            stepAsymData=adaptData.getParamInCond('stepLengthAsym','catch');
            tmcatchData=adaptData.getParamInCond(params,'catch');
            if isempty(tmcatchData)
                newtmcatchData=NaN(1,length(params));
                newStepAsymData=NaN;
            elseif size(tmcatchData,1)<catchNumPts
                newtmcatchData=nanmean(tmcatchData);
                newStepAsymData=nanmean(stepAsymData);
            else
                [newStepAsymData,~]=bin_dataV1(stepAsymData,catchNumPts);
                [newtmcatchData,~]=bin_dataV1(tmcatchData,catchNumPts);
            end        
            [~,maxLoc]=max(abs(newStepAsymData),[],1);
            %ind=sub2ind(size(newtmcatchData),maxLoc*ones(1,length(params)),1:length(params));
            tmcatch=[tmcatch; newtmcatchData(maxLoc,:)];
            stepsymmetryCatch=[stepsymmetryCatch; newStepAsymData(maxLoc,:).*ones(1,size(tmcatch, 2))];
           %%% ~~~~~~~~
         
        %calculate TM steady state #2 
        tmsteady2Data=adaptData.getParamInCond(params,'re-adaptation');
        tmsteady2=[tmsteady2;nanmean(tmsteady2Data((end-5)-steadyNumPts+1:(end-5),:))];
        
        %Magnitude adapted in the first adaptation
        MagAdapt1=[MagAdapt1; nanmean(tmsteady2Data((end-5)-(steadyNumPts)+1:(end-5),:))-nanmean(TMstartData(1:5,:), 1)]
        
        
        %First 5 strides of Transfer
        transferData=adaptData.getParamInCond(params,'OG post');
         
        %%% ~~~~~~~~
        %Band Exclusion
        %Here is where I have the strides that I would have excluded
        %before, need to recalculate the band bass, but its a start.
        %adaptData.getParamInCond({'WhatsUp'},'OG post')
        DStimes = [adaptData.getParamInCond({'doubleSupportSlow'},'OG post'); adaptData.getParamInCond({'doubleSupportFast'},'OG post')];
        SwingTimes = [adaptData.getParamInCond({'swingTimeSlow'},'OG post'); adaptData.getParamInCond({'swingTimeFast'},'OG post')];
        
        if nanmean(DStimes)<0
        DSthresh = -1.5.*nanmean(DStimes);
        Swingthresh = -1.5.*nanmean(SwingTimes);
        else
        DSthresh = 1.5.*nanmean(DStimes);
        Swingthresh = 1.5.*nanmean(SwingTimes);
        end
        
        %DSCorridor=[1.5.*bin_dataV2(DStimes,5)'; .5.*bin_dataV2(DStimes,5)'];
        %SCorridor=[1.5.*bin_dataV2(SwingTimes,5)'; .5.*bin_dataV2(SwingTimes,5)'];
        
        BandData=adaptData.getParamInCond({'doubleSupportSlow', 'doubleSupportFast', 'swingTimeSlow', 'swingTimeFast'},'OG post');
        hipster=adaptData.getParamInCond({'WhatsUp'},'OG post');
        
        BadIndice=[];
%         for lego=1:length(BandData)
%             if BandData(lego, 1)>DSthresh || BandData(lego, 2)>DSthresh ||  BandData(lego, 3)>Swingthresh ||  BandData(lego, 4)>Swingthresh
%                 BadIndice=[BadIndice; lego];
%             end
%         end
%         for lego=1:length(hipster)
%             if hipster(lego)~=0
%                 BadIndice=[BadIndice; lego];
%             end
%         end
% 
%         transferData(BadIndice, :)=[];
         
         if isempty(transferData)==1
             ogafter=[ogafter; nan(1, size(params, 2))];
         else
        ogafter=[ogafter; mean(transferData(1:transientNumPts, :))];
         end
         %%% ~~~~~~~~
         
%         %calculate TM after-effects same as transfer
%         tmafterData=adaptData.getParamInCond(params,'TM post');
%         tmafter=[tmafter; mean(tmafterData(1:transientNumPts, :))];
%         
%         
               %%% ~~~~~~~~
%             calculate TM after-effects as mean value during strides which caused a
%             maximum deviation from zero in step length asymmetry during         
            WstepAsymData=adaptData.getParamInCond('stepLengthAsym','TM post');
            tmafterData=adaptData.getParamInCond(params,'TM post');
            %tmafterData=tmafterData(1:10,:);
            if isempty(tmafterData)
                newtmafterData=NaN(1,length(params));
                newWSAData=NaN;
            elseif size(tmafterData,1)<transientNumPts
                newtmafterData=nanmean(tmafterData);
                newWSAData=nanmean(WstepAsymData);
            else
                [newWSAData,~]=bin_dataV1(WstepAsymData(1:25,:),transientNumPts);
                [newtmafterData,~]=bin_dataV1(tmafterData(1:25,:),transientNumPts);
            end        
            [~,maxLoc]=max(abs(newWSAData),[],1);
            %ind=sub2ind(size(newtmafterData),maxLoc*ones(1,length(params)),1:length(params));
            tmafter=[tmafter; newtmafterData(maxLoc,:)];
          
            %%% ~~~~~~~~
        
        if strcmp(groups{g}, 'OA')==1
            TempTM=[0.133260439601791;0.264262815528979;0.270954600503630;0.232559310763997;0.272167427847111;0.258885239250101;0.371202830537235;0.329968622625190;0.286322713662362;0.155804463282275;0.392077997738293];
            TempTM = (repmat(TempTM,1,size(params, 2)));
        elseif strcmp(groups{g}, 'YA')==1
            TempTM=[[0.279516778950570;0.192298366835098;0.323923392840122;0.358500825965512;0.313055348281848;0.208594489243970;0.438429838146457;0.209597622046027]];
            TempTM = (repmat(TempTM,1,size(params, 2)));
        elseif strcmp(groups{g}, 'OASV')==1
            TempTM=[[0.169267873440147;0.325903711305213;0.322871486313800;0.178460389689225;0.231949816538597;0.195546372835415;0.423510938065962;0.379756841456521]];
            TempTM = (repmat(TempTM,1,size(params, 2)));
        elseif strcmp(groups{g}, 'YASV')==1
            TempTM=[[0.214700576038849;0.241866161697763;0.369100726301672;0.228612381029869;0.441094174658979;0.231270380160492]];
            TempTM = (repmat(TempTM,1,size(params, 2)));
        end
        
        %calculate relative after-effects
% % %         transfer=[transfer; 100*(ogafter./stepsymmetryCatch)];
% % %         washout=[washout; 100*(tmafter./stepsymmetryCatch)];
% % %         %transfer=[transfer; 100*(ogafter./tmcatch)];%Origonal
        
        %Now used
        transfer=[transfer; 100*(ogafter(s,:)./TempTM(s,:))];
        washout=[washout; 100*(tmafter(s,:)./TempTM(s,:))];
        
        transfer2=[transfer2; 100*(ogafter(s,:)./tmsteady2(s,:))];
        washout2=[washout2; 100*(tmafter(s,:)./tmsteady2(s,:))];
        
%         %calculate the average "remembering" during "TM post"
%         epoch=find(ismember(test, 'TM post')==1);
%         tt=adaptData.metaData.trialsInCondition{epoch};
%         rememberDataB1=adaptData.getParamInTrial(params,tt(1));
%         rememberDataB2=adaptData.getParamInTrial(params,tt(2));
%         rememberDataB3=adaptData.getParamInTrial(params,tt(3));
        
%         perR=[(nanmean(rememberDataB1(end-4:end,:))-nanmean(rememberDataB2(1:5,:)))./nanmean(rememberDataB1(end-4:end,:));...
%             (nanmean(forgetB2Data(end-4:end,:))-nanmean(rememberDataB3(1:5,:)))./nanmean(forgetB2Data(end-4:end,:))];
%         remember=[remember; 100*(nanmean(perR))];
%         
        
    end
        results.TMstart.avg(end+1,:)=nanmean(TMstart,1);
    results.TMstart.sd(end+1,:)=nanstd(TMstart,1)/sqrt(length(TMstart));
    results.TMstart.indiv.(groups{g})=TMstart;
    
    results.Forget.avg(end+1,:)=nanmean(forget,1);
    results.Forget.sd(end+1,:)=nanstd(forget,1)/sqrt(length(forget));
    results.Forget.indiv.(groups{g})=forget;
    
    results.AVGForget.avg(end+1,:)=nanmean(AVGforget,1);
    results.AVGForget.sd(end+1,:)=nanstd(AVGforget,1)/sqrt(length(AVGforget));
    results.AVGForget.indiv.(groups{g})=AVGforget;
    
    results.PerForget.avg(end+1,:)=nanmean(perforget,1);
    results.PerForget.sd(end+1,:)=nanstd(perforget,1)/sqrt(length(perforget));
    results.PerForget.indiv.(groups{g})=perforget;
    results.PerForget.indivAll.(groups{g})=perALL;
    
    results.TMsteady1.avg(end+1,:)=nanmean(tmsteady1,1);
    results.TMsteady1.sd(end+1,:)=nanstd(tmsteady1,1)/sqrt(length(tmsteady1));
    results.TMsteady1.indiv.(groups{g})=tmsteady1;
    
    results.catch.avg(end+1,:)=nanmean(tmcatch,1);
    results.catch.sd(end+1,:)=nanstd(tmcatch,1)/sqrt(length(tmsteady1));
    results.catch.indiv.(groups{g})=tmcatch;
    
    results.TMsteady2.avg(end+1,:)=nanmean(tmsteady2,1);
    results.TMsteady2.sd(end+1,:)=nanstd(tmsteady2,1)/sqrt(length(tmsteady1));
    results.TMsteady2.indiv.(groups{g})=tmsteady2;
    
    results.TMafter.avg(end+1,:)=nanmean(tmafter,1);
    results.TMafter.sd(end+1,:)=nanstd(tmafter,1)/sqrt(length(tmsteady1));
    results.TMafter.indiv.(groups{g})=tmafter;    
    
    results.OGafter.avg(end+1,:)=nanmean(ogafter,1);
    results.OGafter.sd(end+1,:)=nanstd(ogafter,1)/sqrt(length(tmsteady1));
    results.OGafter.indiv.(groups{g})=ogafter;
    
    results.Transfer.avg(end+1,:)=nanmean(transfer,1);
    results.Transfer.sd(end+1,:)=nanstd(transfer,1)/sqrt(length(tmsteady1));
    results.Transfer.indiv.(groups{g})=transfer;
    
    results.Washout.avg(end+1,:)=nanmean(washout,1);
    results.Washout.sd(end+1,:)=nanstd(washout,1)/sqrt(length(tmsteady1));
    results.Washout.indiv.(groups{g})=washout;
    
    results.Transfer2.avg(end+1,:)=nanmean(transfer2,1);
    results.Transfer2.sd(end+1,:)=nanstd(transfer2,1)/sqrt(length(tmsteady1));
    results.Transfer2.indiv.(groups{g})=transfer2;
    
    results.Washout2.avg(end+1,:)=nanmean(washout2,1);
    results.Washout2.sd(end+1,:)=nanstd(washout2,1)/sqrt(length(tmsteady1));
    results.Washout2.indiv.(groups{g})=washout2;
    
    results.Remember.avg(end+1,:)=nanmean(remember,1);
    results.Remember.sd(end+1,:)=nanstd(remember,1)/sqrt(length(tmsteady1));
    results.Remember.indiv.(groups{g})=remember;
    

    
    results.MagAdapt1.avg(end+1,:)=nanmean(MagAdapt1,1);
    results.MagAdapt1.sd(end+1,:)=nanstd(MagAdapt1,1)/sqrt(length(MagAdapt1));
    results.MagAdapt1.indiv.(groups{g})=MagAdapt1;
end

%plot stuff
epochs=fields(results);

% %plot first five epochs
% numPlots=5*length(params); 
% ah=optimizedSubPlot(numPlots,length(params),5,'ltr');
% i=1;
% for p=1:length(params)
%     limy=[];
%     for t=[1:5]   
%         axes(ah(i))
%         hold on        
%         for b=1:ngroups
%             bar(b,results.(epochs{t}).avg(b,p),'facecolor',ColorOrder(b,:));
%             if nargin>3 && ~isempty(indivFlag)
%                 plot(b,results.(epochs{t}).indiv.(groups{b})(:,p),'k*')
%             end
%         end
%         errorbar(results.(epochs{t}).avg(:,p),results.(epochs{t}).sd(:,p),'.','LineWidth',2,'Color','k')
%         set(gca,'Xtick',1:ngroups,'XTickLabel',groups,'fontSize',12)
%         axis tight
%         limy=[limy get(gca,'Ylim')];
%         ylabel(params{p})
%         title(epochs{t})
%         i=i+1;
%         
%     end
%     set(ah(p*5-4:p*5),'Ylim',[min(limy) max(limy)])
%     %set(ah(p*5-4:p*5),'Ylim',[0 200])
% end


% %plot last four epochs
% numPlots=5*length(params);
% ah=optimizedSubPlot(numPlots,length(params),5,'ltr');
% i=1;
% for p=1:length(params)
%     limy=[];
%     for t=8:12
%         axes(ah(i))
%         hold on        
%         for b=1:ngroups
%             bar(b,results.(epochs{t}).avg(b,p),'facecolor',ColorOrder(b,:));            
%         end
%         errorbar(results.(epochs{t}).avg(:,p),results.(epochs{t}).sd(:,p),'.','LineWidth',2,'Color','k')
%         set(gca,'Xtick',1:ngroups,'XTickLabel',groups,'fontSize',12)
%         axis tight
%         limy=[limy get(gca,'Ylim')];
%         ylabel(params{p})
%         title(epochs{t})
%         i=i+1;
%     end
%     set(ah(p*5-3:p*5),'Ylim',[min(limy) max(limy)])
% end

%plot first five epochs
numPlots=6*length(params); 
ah=optimizedSubPlot(numPlots,length(params),6,'ltr');
set(ah,'defaultTextFontName', 'Arial')
i=1;limymin=[];limymax=[];
for p=1:length(params)
    
    for t=[1:6]   
        axes(ah(i))
        hold on        
%         for b=1:ngroups
%             %something(b) = subplot(b);
%             %bb=gca
%             bar(b,results.(epochs{t}).avg(b,p), 'facecolor',ColorOrder(b,:));
%             %bar(b,results.(epochs{t}).avg(b,p))%,'facecolor',ColorOrder(b,:));
%             if nargin>3 && ~isempty(indivFlag)
%                 plot(b,results.(epochs{t}).indiv.(groups{b})(:,p),'k*')
%             end
%         end
%         errorbar(results.(epochs{t}).avg(:,p),results.(epochs{t}).sd(:,p),'.','LineWidth',2,'Color','k')
%         set(gca,'Xtick',1:ngroups,'XTickLabel',groups,'fontSize',12)

line([.5 4.5], [0 0], 'Color','k')
      for b=1:ngroups
       
        bar(b,results.(epochs{t}).avg(b,p),'facecolor',ColorOrder(b,:));
            %plot(b.*ones(1, size(results.PerForget.indiv.(groups{b}), 1)), results.PerForget.indiv.(groups{b})(:,p), 'k.')
        errorbar(b, results.(epochs{t}).avg(b,p),results.(epochs{t}).sd(b,p),'.','LineWidth',2,'Color','k')

%%% For Plotting only the Color Coded Error Bars 
%         errorbar(b, results.(epochs{t}).avg(b,p),results.(epochs{t}).sd(b,p),'.','LineWidth',2,'Color',ColorOrder(b,:))
       
       
      end
% line([1 3], [results.(epochs{t}).avg(1,p) results.(epochs{t}).avg(3,p)], 'Marker','.','LineStyle','--', 'Color','k')
% line([2 4], [results.(epochs{t}).avg(2,p) results.(epochs{t}).avg(4,p)], 'Marker','.','LineStyle',':', 'Color','k')

%set(gca,'Xtick',1:ngroups,'XTickLabel',[1 2 1 2],'fontSize',12)
        set(gca,'Xtick',[1.5 3.5],'XTickLabel',[{'Old'} {'Young'}],'fontSize',12)

        axis tight
        %limy=[limy get(gca,'Ylim')];
        temp=[get(gca,'Ylim')];
        limymin(p, t)=temp(1);
        limymax(p, t)=temp(2);
        ylabel(params{p})
        title(epochs{t})
        i=i+1;
        
    end
    %set(ah(p*5-4:p*5),'Ylim',[min(limy) max(limy)])
    %set(ah(p*5-4:p*5),'Ylim',[0 200])
    %set(ah(p*5-4:p*5),'Xlim',[0.5 4.5])
end
set(ah(1:6*p),'Xlim',[0.5 4.5])

set(ah(1),'Ylim',[min([limymin(1,1) limymin(2,1) limymin(3,1) limymin(4,1)]) max([limymax(1,1) limymax(2,1) limymax(3,1) limymax(4,1)])])
set(ah(7),'Ylim',[min([limymin(1,1) limymin(2,1) limymin(3,1) limymin(4,1)]) max([limymax(1,1) limymax(2,1) limymax(3,1) limymax(4,1)])])
set(ah(13),'Ylim',[min([limymin(1,1) limymin(2,1) limymin(3,1) limymin(4,1)]) max([limymax(1,1) limymax(2,1) limymax(3,1) limymax(4,1)])])
set(ah(19),'Ylim',[min([limymin(1,1) limymin(2,1) limymin(3,1) limymin(4,1)]) max([limymax(1,1) limymax(2,1) limymax(3,1) limymax(4,1)])])

set(ah(2),'Ylim',[min([limymin(1,2) limymin(2,2) limymin(3,2) limymin(4,2)]) max([limymax(1,2) limymax(2,2) limymax(3,2) limymax(4,2)])])
set(ah(8),'Ylim',[min([limymin(1,2) limymin(2,2) limymin(3,2) limymin(4,2)]) max([limymax(1,2) limymax(2,2) limymax(3,2) limymax(4,2)])])
set(ah(14),'Ylim',[min([limymin(1,2) limymin(2,2) limymin(3,2) limymin(4,2)]) max([limymax(1,2) limymax(2,2) limymax(3,2) limymax(4,2)])])
set(ah(20),'Ylim',[min([limymin(1,2) limymin(2,2) limymin(3,2) limymin(4,2)]) max([limymax(1,2) limymax(2,2) limymax(3,2) limymax(4,2)])])

set(ah(3),'Ylim',[min([limymin(1,3) limymin(2,3) limymin(3,3) limymin(4,3)]) max([limymax(1,3) limymax(2,3) limymax(3,3) limymax(4,3)])])
set(ah(9),'Ylim',[min([limymin(1,3) limymin(2,3) limymin(3,3) limymin(4,3)]) max([limymax(1,3) limymax(2,3) limymax(3,3) limymax(4,3)])])
set(ah(15),'Ylim',[min([limymin(1,3) limymin(2,3) limymin(3,3) limymin(4,3)]) max([limymax(1,3) limymax(2,3) limymax(3,3) limymax(4,3)])])
set(ah(21),'Ylim',[min([limymin(1,3) limymin(2,3) limymin(3,3) limymin(4,3)]) max([limymax(1,3) limymax(2,3) limymax(3,3) limymax(4,3)])])

    
% 
set(ah(5),'Ylim',[-.1 .3])
set(ah(11),'Ylim',[-.01 .11])
set(ah(17),'Ylim',[-.3 .02])
set(ah(23),'Ylim',[-.3 .02])

% set(ah(4),'Ylim',[-.1 .3])
% set(ah(10),'Ylim',[-.01 .11])
% set(ah(16),'Ylim',[-.3 .02])
% set(ah(22),'Ylim',[-.3 .02])

% set(ah(5),'Ylim',[-.08537 0.2257])
% set(ah(11),'Ylim',[-.01433 0.1097])
% set(ah(17),'Ylim',[-.2983 0.007779])
% set(ah(23),'Ylim',[-.2935 0.02669])
% 
% set(ah(4),'Ylim',[-.08537 0.2257])
% set(ah(10),'Ylim',[-.01433 0.1097])
% set(ah(16),'Ylim',[-.2983 0.007779])
% set(ah(22),'Ylim',[-.2935 0.02669])

% set(ah(5),'Ylim',[min([limymin(1,5) limymin(1,6)]) max([limymax(1,5) limymax(1,6)])])
% set(ah(11),'Ylim',[min([limymin(2,5) limymin(2,6)]) max([limymax(2,5) limymax(2,6)])])
% set(ah(17),'Ylim',[min([ limymin(3,5) limymin(3,6)]) max([limymax(3,5) limymax(3,6)])])
% set(ah(23),'Ylim',[min([limymin(4,5) limymin(4,6)]) max([limymax(4,5) limymax(4,6)])])

% set(ah(6),'Ylim',[min([limymin(1,5) limymin(1,6)]) max([limymax(1,5) limymax(1,6)])])
% set(ah(12),'Ylim',[min([limymin(2,5) limymin(2,6)]) max([limymax(2,5) limymax(2,6)])])
% set(ah(18),'Ylim',[min([limymin(3,5) limymin(3,6)]) max([limymax(3,5) limymax(3,6)])])
% set(ah(24),'Ylim',[min([limymin(4,5) limymin(4,6)]) max([limymax(4,5) limymax(4,6)])])


set(ah(4),'Ylim',[min([limymin(1,4) limymin(2,4) limymin(3,4) limymin(4,4)]) max([limymax(1,4) limymax(2,4) limymax(3,4) limymax(4,4)])])
set(ah(10),'Ylim',[min([limymin(1,4) limymin(2,4) limymin(3,4) limymin(4,4)]) max([limymax(1,4) limymax(2,4) limymax(3,4) limymax(4,4)])])
set(ah(16),'Ylim',[min([limymin(1,4) limymin(2,4) limymin(3,4) limymin(4,4)]) max([limymax(1,4) limymax(2,4) limymax(3,4) limymax(4,4)])])
set(ah(22),'Ylim',[min([limymin(1,4) limymin(2,4) limymin(3,4) limymin(4,4)]) max([limymax(1,4) limymax(2,4) limymax(3,4) limymax(4,4)])])
% 
% set(ah(5),'Ylim',[min([limymin(1,5) limymin(2,5) limymin(3,5) limymin(4,5)]) max([limymax(1,5) limymax(2,5) limymax(3,5) limymax(4,5)])])
% set(ah(11),'Ylim',[min([limymin(1,5) limymin(2,5) limymin(3,5) limymin(4,5)]) max([limymax(1,5) limymax(2,5) limymax(3,5) limymax(4,5)])])
% set(ah(17),'Ylim',[min([limymin(1,5) limymin(2,5) limymin(3,5) limymin(4,5)]) max([limymax(1,5) limymax(2,5) limymax(3,5) limymax(4,5)])])
% set(ah(23),'Ylim',[min([limymin(1,5) limymin(2,5) limymin(3,5) limymin(4,5)]) max([limymax(1,5) limymax(2,5) limymax(3,5) limymax(4,5)])])

set(ah(6),'Ylim',[min([limymin(1,6) limymin(2,6) limymin(3,6) limymin(4,6)]) max([limymax(1,6) limymax(2,6) limymax(3,6) limymax(4,6)])])
set(ah(12),'Ylim',[min([limymin(1,6) limymin(2,6) limymin(3,6) limymin(4,6)]) max([limymax(1,6) limymax(2,6) limymax(3,6) limymax(4,6)])])
set(ah(18),'Ylim',[min([limymin(1,6) limymin(2,6) limymin(3,6) limymin(4,6)]) max([limymax(1,6) limymax(2,6) limymax(3,6) limymax(4,6)])])
set(ah(24),'Ylim',[min([limymin(1,6) limymin(2,6) limymin(3,6) limymin(4,6)]) max([limymax(1,6) limymax(2,6) limymax(3,6) limymax(4,6)])])

limymin=[]; limymax=[];

%plot last four epochs
numPlots=6*length(params);
ah=optimizedSubPlot(numPlots,length(params),6,'ltr');
i=1;limymin=[];limymax=[];
for p=1:length(params)
    
% % %     for t=[7:12]
% % %         axes(ah(i))
% % %         hold on        
% % % %         for b=1:ngroups
% % % %             bar(b,results.(epochs{t}).avg(b,p),'facecolor',ColorOrder(b,:));            
% % % %         end
% % % %         errorbar(results.(epochs{t}).avg(:,p),results.(epochs{t}).sd(:,p),'.','LineWidth',2,'Color','k')
% % % %         set(gca,'Xtick',1:ngroups,'XTickLabel',groups,'fontSize',12)
% % % line([.5 4.5], [0 0], 'Color','k')
% % %         for b=1:ngroups
% % %         %%%errorbar(results.(epochs{t}).avg(b,p),results.(epochs{t}).sd(b,p),'.','LineWidth',2,'Color',ColorOrder(b,:))
% % %         errorbar(b, results.(epochs{t}).avg(b,p),results.(epochs{t}).sd(b,p),'.','LineWidth',2,'Color',ColorOrder(b,:))
% % %        end
% % % line([1 3], [results.(epochs{t}).avg(1,p) results.(epochs{t}).avg(3,p)], 'Marker','.','LineStyle','--', 'Color','k')
% % % line([2 4], [results.(epochs{t}).avg(2,p) results.(epochs{t}).avg(4,p)], 'Marker','.','LineStyle',':', 'Color','k')
% % 
% % %set(gca,'Xtick',1:ngroups,'XTickLabel',[1 2 1 2],'fontSize',12)

    for t=[7:12]
        axes(ah(i))
        hold on  

        for b=1:ngroups
bar(b,results.(epochs{t}).avg(b,p),'facecolor',ColorOrder(b,:));

        errorbar(b, results.(epochs{t}).avg(b,p),results.(epochs{t}).sd(b,p),'.','LineWidth',2,'Color','k')
    end
        set(gca,'Xtick',[1.5 3.5],'XTickLabel',[{'Old'} {'Young'}],'fontSize',12)
        
        axis tight
        temp=[get(gca,'Ylim')];
        limymin(p, t-6)=temp(1);
        limymax(p, t-6)=temp(2);
        ylabel(params{p})
        title(epochs{t})
        i=i+1;
    end
    %set(ah(p*5-3:p*5),'Ylim',[min(limy) max(limy)])
    
end
set(ah(1:6*p),'Xlim',[0.5 4.5])
display('fin')

set(ah(1),'Ylim',[-.03 .09])
set(ah(7),'Ylim',[-.03 .09])
set(ah(13),'Ylim',[-.03 .09])
set(ah(19),'Ylim',[-.03 .09])


% set(ah(1),'Ylim',[min([limymin(1,1) limymin(2,1) limymin(3,1) limymin(4,1)]) max([limymax(1,1) limymax(2,1) limymax(3,1) limymax(4,1)])])
% set(ah(7),'Ylim',[min([limymin(1,1) limymin(2,1) limymin(3,1) limymin(4,1)]) max([limymax(1,1) limymax(2,1) limymax(3,1) limymax(4,1)])])
% set(ah(13),'Ylim',[min([limymin(1,1) limymin(2,1) limymin(3,1) limymin(4,1)]) max([limymax(1,1) limymax(2,1) limymax(3,1) limymax(4,1)])])
% set(ah(19),'Ylim',[min([limymin(1,1) limymin(2,1) limymin(3,1) limymin(4,1)]) max([limymax(1,1) limymax(2,1) limymax(3,1) limymax(4,1)])])

set(ah(2),'Ylim',[-.03 .25])
set(ah(8),'Ylim',[-.03 .25])
set(ah(14),'Ylim',[-.03 .25])
set(ah(20),'Ylim',[-.03 .25])

% set(ah(2),'Ylim',[min([limymin(1,2) limymin(2,2) limymin(3,2) limymin(4,2)]) max([limymax(1,2) limymax(2,2) limymax(3,2) limymax(4,2)])])
% set(ah(8),'Ylim',[min([limymin(1,2) limymin(2,2) limymin(3,2) limymin(4,2)]) max([limymax(1,2) limymax(2,2) limymax(3,2) limymax(4,2)])])
% set(ah(14),'Ylim',[min([limymin(1,2) limymin(2,2) limymin(3,2) limymin(4,2)]) max([limymax(1,2) limymax(2,2) limymax(3,2) limymax(4,2)])])
% set(ah(20),'Ylim',[min([limymin(1,2) limymin(2,2) limymin(3,2) limymin(4,2)]) max([limymax(1,2) limymax(2,2) limymax(3,2) limymax(4,2)])])

set(ah(3),'Ylim',[min([limymin(1,3) limymin(2,3) limymin(3,3) limymin(4,3)]) max([limymax(1,3) limymax(2,3) limymax(3,3) limymax(4,3)])])
set(ah(9),'Ylim',[min([limymin(1,3) limymin(2,3) limymin(3,3) limymin(4,3)]) max([limymax(1,3) limymax(2,3) limymax(3,3) limymax(4,3)])])
set(ah(15),'Ylim',[min([limymin(1,3) limymin(2,3) limymin(3,3) limymin(4,3)]) max([limymax(1,3) limymax(2,3) limymax(3,3) limymax(4,3)])])
set(ah(21),'Ylim',[min([limymin(1,3) limymin(2,3) limymin(3,3) limymin(4,3)]) max([limymax(1,3) limymax(2,3) limymax(3,3) limymax(4,3)])])

set(ah(4),'Ylim',[min([limymin(1,4) limymin(2,4) limymin(3,4) limymin(4,4)]) max([limymax(1,4) limymax(2,4) limymax(3,4) limymax(4,4)])])
set(ah(10),'Ylim',[min([limymin(1,4) limymin(2,4) limymin(3,4) limymin(4,4)]) max([limymax(1,4) limymax(2,4) limymax(3,4) limymax(4,4)])])
set(ah(16),'Ylim',[min([limymin(1,4) limymin(2,4) limymin(3,4) limymin(4,4)]) max([limymax(1,4) limymax(2,4) limymax(3,4) limymax(4,4)])])
set(ah(22),'Ylim',[min([limymin(1,4) limymin(2,4) limymin(3,4) limymin(4,4)]) max([limymax(1,4) limymax(2,4) limymax(3,4) limymax(4,4)])])

set(ah(5),'Ylim',[min([limymin(1,5) limymin(2,5) limymin(3,5) limymin(4,5)]) max([limymax(1,5) limymax(2,5) limymax(3,5) limymax(4,5)])])
set(ah(11),'Ylim',[min([limymin(1,5) limymin(2,5) limymin(3,5) limymin(4,5)]) max([limymax(1,5) limymax(2,5) limymax(3,5) limymax(4,5)])])
set(ah(17),'Ylim',[min([limymin(1,5) limymin(2,5) limymin(3,5) limymin(4,5)]) max([limymax(1,5) limymax(2,5) limymax(3,5) limymax(4,5)])])
set(ah(23),'Ylim',[min([limymin(1,5) limymin(2,5) limymin(3,5) limymin(4,5)]) max([limymax(1,5) limymax(2,5) limymax(3,5) limymax(4,5)])])

set(ah(6),'Ylim',[min([limymin(1,6) limymin(2,6) limymin(3,6) limymin(4,6)]) max([limymax(1,6) limymax(2,6) limymin(3,6) limymin(4,6)])])
set(ah(12),'Ylim',[min([limymin(1,6) limymin(2,6) limymin(3,6) limymin(4,6)]) max([limymax(1,6) limymax(2,6) limymax(3,6) limymin(4,6)])])
set(ah(18),'Ylim',[min([limymin(1,6) limymin(2,6) limymin(3,6) limymin(4,6)]) max([limymax(1,6) limymax(2,6) limymax(3,6) limymin(4,6)])])
set(ah(24),'Ylim',[min([limymin(1,6) limymin(2,6) limymin(3,6) limymin(4,6)]) max([limymax(1,6) limymax(2,6) limymax(3,6) limymin(4,6)])])