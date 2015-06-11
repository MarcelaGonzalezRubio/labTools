function []=BFerror2rev2(subject) 

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
adaptData = adaptData.removeBias;
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
XFast=adaptData.getParamInCond('XFast',condition(i)); %fAnkPos(indSHS)
betaFast=adaptData.getParamInCond('betaFast',condition(i)); %FTO
alphaFast=adaptData.getParamInCond('alphaFast',condition(i)); %FHS
RatioFastPos=adaptData.getParamInCond('RFastPos',condition(i)); %FTO/FHS
RFastPosSHS=adaptData.getParamInCond('RFastPosSHS',condition(i));%fAnkPos(SHS)/FHS

betaSlow=adaptData.getParamInCond('betaSlow',condition(i)); %STO
alphaSlow=adaptData.getParamInCond('alphaTemp',condition(i)); %SHS
RatioSlowPos=adaptData.getParamInCond('RSlowPos',condition(i)); %STO/SHS
XSlow=adaptData.getParamInCond('XSlow',condition(i)); %sAnkPos(indFHS)
RSlowPosFHS=adaptData.getParamInCond('RSlowPosFHS',condition(i)); %sAnkPos(FHS)/SHS

if strcmp(condition{i},'Gradual adaptation') || strcmp(condition{i},'Gradual Adaptation') || strcmp(condition{i},'adaptation')
RatioFastPos=[];
RatioSlowPos=[];
RatioFastFHS=[];
RatioSlowSHS=[];
RatioFastPos=adaptData.getParamInCond('RSlowPos','TM base');
RatioSlowPos=adaptData.getParamInCond('RFastPos','TM base');
RSlowPosFHS=adaptData.getParamInCond('RSlowPosFHS','TM base');
RFastPosSHS=adaptData.getParamInCond('RFastPosSHS','TM base');
end

meanFastPos=mean(RatioFastPos);
meanSlowPos=mean(RatioSlowPos);
meanFastSHS=mean(RFastPosSHS);
meanSlowFHS=mean(RSlowPosFHS);

%Error calculation e=|beta|-|R*alpha|
for n=1:1:length(betaFast)
    EFastbeta(n,1)=abs(betaFast(n,1))-abs(meanFastPos*alphaFast(n,1));
    ESlowbeta(n,1)=abs(betaSlow(n,1))-abs(meanSlowPos*alphaSlow(n,1)); 
    
    EFastX(n,1)=abs(XFast(n,1))-abs(meanFastSHS*alphaFast(n,1));
    ESlowX(n,1)=abs(XSlow(n,1))-abs(meanSlowFHS*alphaSlow(n,1));
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
