%Simulation of feedback 
function []=Feedback_test(subject)

load([subject 'params.mat'])
load([subject 'RAW.mat'])

condition= adaptData.metaData.conditionName;
% This is new ~~~~~~~~~~~~~~~~~
     %Needed in case subject did not perform one of the conditions
        %in the condition list
condition=condition(find(~cellfun(@isempty,adaptData.metaData.trialsInCondition)));

alphaFastN=[];
alphaFastNHipmean=[];
alphaFastNWOHip=[]
% alphaFastO=[];
alphaSlowN=[];
alphaSlowNHipmean=[]
alphaSlowNWOHip=[];
% alphaSlowO=[];

for i=1:1:length(condition)

    alphaSlow=adaptData.getParamInCond('alphaTemp',condition(i));
    alphaSlowHipmean=adaptData.getParamInCond('alphaTempHipmean',condition(i));
    alphaSlowWOHip=adaptData.getParamInCond('alphaSlowWOHip',condition(i));
    
    
    alphaFast=adaptData.getParamInCond('alphaFast',condition(i));
    alphaFastHipmean=adaptData.getParamInCond('alphaFastHipmean',condition(i));
    alphaFastWOHip=adaptData.getParamInCond('alphaFastWOHip',condition(i));
    
    
    betaSlow=adaptData.getParamInCond('betaSlow',condition(i));
    betaSlowHipmean=adaptData.getParamInCond('betaSlowHipmean',condition(i));
    betaSlowWOHip=adaptData.getParamInCond('betaSlowWOHip',condition(i));
    
    betaFast=adaptData.getParamInCond('betaFast',condition(i));
    betaFastHipmean=adaptData.getParamInCond('betaFastHipmean',condition(i));
    betaFastWOHip=adaptData.getParamInCond('betaFastWOHip',condition(i));
    
    xSlow=adaptData.getParamInCond('xSlow',condition(i));
    xSlowHipmean=adaptData.getParamInCond('xSlowHipmean',condition(i));
    xSlowWOHip=adaptData.getParamInCond('xSlowWOHip',condition(i));
    
    xFast=adaptData.getParamInCond('xFast',condition(i));
    xFastHipmean=adaptData.getParamInCond('xFastHipmean',condition(i));
    xFastWOHip=adaptData.getParamInCond('xFastWOHip',condition(i));

alphaFastNew=alphaFast;
alphaFastNewHipmean=alphaFastHipmean;
alphaFastNewWOHip=alphaFastWOHip;

alphaSlowNew=alphaSlow;
alphaSlowNewHipmean=alphaSlowHipmean;
alphaSlowNewWOHip=alphaSlowWOHip;

% alphaFastOld=alphaFast;
% alphaSlowOld=alphaSlow;


% if strcmp(condition{i},'TM base')
%     RSlow=adaptData.getParamInCond('RSlowPos',condition(i));
%     RFast=adaptData.getParamInCond('RFastPos',condition(i));
%     RSlow=nanmean(RSlow);
%     RFast=nanmean(RFast);
%     
% end


if strcmp(condition{i},'Gradual adaptation') || strcmp(condition{i},'Re-adaptation') || strcmp(condition{i},'adaptation')
    
%     alphaFastNew=alphaFast;
%     alphaSlowNew=alphaSlow;
    
    for n=2:1:length(alphaSlow)
        
        alphaSlowNew(n)=alphaFastNew(n-1)-xSlow(n-1)+xFast(n-1);
        alphaSlowNewHipmean(n)=alphaFastNewHipmean(n-1)-xSlowHipmean(n-1)+xFastHipmean(n-1);
        alphaSlowNewWOHip(n)=alphaFastNewWOHip(n-1)-xSlowWOHip(n-1)+xFastWOHip(n-1);
        
        alphaFastNew(n+1)=alphaSlowNew(n)+betaSlow(n)-betaFast(n);
        alphaFastNewHipmean(n+1)=alphaSlowNewHipmean(n)+betaSlowHipmean(n)-betaFastHipmean(n);
        alphaFastNewWOHip(n+1)=alphaSlowNewWOHip(n)+betaSlowWOHip(n)-betaFastWOHip(n);
        
       
    end
    
%     alphaFastOld=alphaFast;
%     alphaSlowOld=alphaSlow;
    
%     n=0;
%     
%     for n=2:1:length(alphaSlow)
%         alphaFastOld(n)=alphaFastOld(n-1)+(1/(1+RFast))*(abs(betaFast(n-1))-RFast*(alphaFastOld(n-1)));
%         alphaSlowOld(n)=alphaSlowOld(n-1)+(1/(1+RSlow))*(abs(betaSlow(n-1))-RSlow*(alphaSlowOld(n-1)));
%     end
    alphaFastNew=alphaFastNew(1:end-1);
    alphaFastNewHipmean=alphaFastNewHipmean(1:end-1);
    alphaFastNewWOHip=alphaFastNewWOHip(1:end-1); 
end

alphaFastN=[alphaFastN; alphaFastNew];
alphaFastNHipmean=[alphaFastNHipmean; alphaFastNewHipmean];
alphaFastNWOHip=[alphaFastNWOHip; alphaFastNewWOHip];

%alphaFastO=[alphaFastO; alphaFastOld];

alphaSlowN=[alphaSlowN; alphaSlowNew];
alphaSlowNHipmean=[alphaSlowNHipmean; alphaSlowNewHipmean];
alphaSlowNWOHip=[alphaSlowNWOHip; alphaSlowNewWOHip];
%alphaSlowO=[alphaSlowO; alphaSlowOld];

for t=1:1:length(alphaSlowN)
    
    alphaSlowPlusbetaSlow(t)= (alphaSlowN(t)+betaSlow(t))/abs(alphaSlowN(t)-betaSlow(t));
    alphaSlowPlusbetaSlowmean(t)= (alphaSlowNHipmean(t)+betaSlowHipmean(t))/abs(alphaSlowNHipmean(t)-betaSlowHipmean(t));
    alphaSlowPlusbetaSlowWOHip(t)=(((betaSlowWOHip(t)-alphaSlowNWOHip(t))/2)+alphaSlowNWOHip(t))/abs(betaFastWOHip(t));
    
    alphaFastPlusbetaFast(t)=(alphaFastN(t)+betaFast(t))/abs(alphaFastN(t)-betaFast(t));
    alphaFastPlusbetaFastmean(t)=(alphaFastNHipmean(t)+betaFastHipmean(t))/abs(alphaFastNHipmean(t)-betaFastHipmean(t));
    alphaFastPlusbetaFastWOHip(t)=(((betaFastWOHip(t)-alphaFastNWOHip(t))/2)+alphaFastNWOHip(t))/abs(betaFastWOHip(t));
    
    alphaSlowPlusXSlow(t)=(alphaSlowN(t)+xSlow(t))/abs(alphaSlowN(t)-xSlow(t));
    alphaSlowPlusXSlowmean(t)=(alphaSlowNHipmean(t)+xSlowHipmean(t))/abs(alphaSlowNHipmean(t)-betaSlowHipmean(t));
    alphaSlowPlusXSlowWOHip(t)=(((xSlowWOHip(t)-alphaSlowNWOHip(t))/2)+alphaSlowNWOHip(t))/abs(betaFastWOHip(t));
    
    alphaFastPlusXFast(t)=(alphaFastN(t)+xFast(t))/abs(alphaFastN(t)-xFast(t));
    alphaFastPlusXFastmean(t)=(alphaFastNHipmean(t)+xFastHipmean(t))/abs(alphaFastNHipmean(t)-betaFastHipmean(t));
    alphaFastPlusXFastWOHip(t)=(((xFastWOHip(t)-alphaFastNWOHip(t))/2)+alphaFastNWOHip(t))/abs(betaFastWOHip(t));
    
end

end


this=paramData([adaptData.data.Data,alphaFastN,alphaSlowN],[adaptData.data.labels 'alphaFastN' 'alphaSlowN'],adaptData.data.indsInTrial,adaptData.data.trialTypes);
test=adaptationData(rawExpData.metaData,rawExpData.subData,this);
saveloc=[];
save([saveloc subject 'test.mat'],'test'); 