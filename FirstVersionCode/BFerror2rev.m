function []=BFerror2rev(subject) 

load([subject 'params.mat'])
load([subject 'RAW.mat'])
condition= adaptData.metaData.conditionName;
% This is new ~~~~~~~~~~~~~~~~~
     %Needed in case subject did not perform one of the conditions
        %in the condition list
        
   condition=condition(find(~cellfun(@isempty,adaptData.metaData.trialsInCondition)));

ErrorFast1=[];
ErrorSlow1=[];
ErrorFast2=[];
ErrorSlow2=[];

for i=1:1:length(condition)

EFastbeta=[];
ESlowbeta=[];
EFastX=[];
ESlowX=[];
betaFast=[];
alphaFast=[];
RatioFastPos=[];
betaSlow=[];
alphaSlow=[];
RatioSlowPos=[];


%Fast leg   
betaFast=adaptData.getParamInCond('betaFast',condition(i));
alphaFast=adaptData.getParamInCond('alphaFast',condition(i));
RatioFastPos=adaptData.getParamInCond('RFastPos',condition(i));
XFast=adaptData.getParamInCond('XFast',condition(i));
RFastPosFHS=adaptData.getParamInCond('RFastPosFHS',condition(i));
betaSlow=adaptData.getParamInCond('betaSlow',condition(i));
alphaSlow=adaptData.getParamInCond('alphaTemp',condition(i));
RatioSlowPos=adaptData.getParamInCond('RSlowPos',condition(i));
XSlow=adaptData.getParamInCond('XSlow',condition(i));
RSlowPosSHS=adaptData.getParamInCond('RSlowPosSHS',condition(i));

% if  length(condition)==9 & i == 4 
%     
%      RatioFast=adaptData.getParamInCond('RFastPos','TM base W/BF');
% else
%     if length(condition)==7 & i==3
%         RatioFast=RatioFast([1:290],1);
%     else
%         RatioFast=RatioFast;
%     end 
% end

%Slow Leg

% if  length(condition)==9 & i == 4 
%     
%      RatioSlow=adaptData.getParamInCond('RSlowPos','TM base W/BF');
% else
%     if length(condition)==7 & i==3
%         RatioSlow=RatioSlow([1:290],1);
%     else
%         RatioSlow=RatioSlow;
%     end 
% end 
% 


if strcmp(condition{i},'Gradual adaptation') || strcmp(condition{i},'Gradual Adaptation') || strcmp(condition{i},'adaptation')
RatioFastPos=[];
RatioSlowPos=[];
RatioFastFHS=[]
RatioSlowSHS=[];
RatioFastPos=adaptData.getParamInCond('RSlowPos','TM base');
RatioSlowPos=adaptData.getParamInCond('RFastPos','TM base');
RSlowPosSHS=adaptData.getParamInCond('RSlowPosSHS','TM base');
RFastPosFHS=adaptData.getParamInCond('RFastPosFHS','TM base');
end

meanFastPos=mean(RatioFastPos);
meanSlowPos=mean(RatioSlowPos);
meanFastFHS=mean(RFastPosFHS);
meanSlowSHS=mean(RSlowPosSHS);

%Error calculation e=|beta|-|R*alpha|
for n=1:1:length(betaFast)
    EFastbeta(n,1)=abs(betaFast(n,1))-abs(meanFastPos*alphaFast(n,1));
    ESlowbeta(n,1)=abs(betaSlow(n,1))-abs(meanSlowPos*alphaSlow(n,1)); 
    
    EFastX(n,1)=abs(XSlow(n,1))-abs(meanFastFHS*alphaFast(n,1));
    ESlowX(n,1)=abs(XFast(n,1))-abs(meanSlowSHS*alphaSlow(n,1));
end

ErrorFast1=[ErrorFast1;EFastbeta];  
ErrorSlow1=[ErrorSlow1;ESlowbeta];
ErrorFast2=[ErrorFast2;EFastX];
ErrorSlow2=[ErrorSlow2;ESlowX];

BFerror1=ErrorFast1-ErrorSlow1;
BFerror2=ErrorFast2-ErrorSlow2;
end

this=paramData([adaptData.data.Data,ErrorFast1,ErrorSlow1,BFerror1,ErrorFast2,ErrorSlow2,BFerror2],[adaptData.data.labels 'Error1Fast' 'Error1Slow' 'TotalError1' 'Error2Fast' 'Error2Slow' 'TotalError2'],adaptData.data.indsInTrial,adaptData.data.trialTypes);
adaptData=adaptationData(rawExpData.metaData,rawExpData.subData,this);
saveloc=[];
save([saveloc subject 'params.mat'],'adaptData'); 

end