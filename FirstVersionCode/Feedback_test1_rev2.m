%Simulation of feedback 
function []=Feedback_test1_rev2(subject)

load([subject 'params.mat'])
load([subject 'RAW.mat'])

condition= adaptData.metaData.conditionName;
% This is new ~~~~~~~~~~~~~~~~~
     %Needed in case subject did not perform one of the conditions
        %in the condition list
condition=condition(find(~cellfun(@isempty,adaptData.metaData.trialsInCondition)));

alphaFastNew=[];
alphaFastOldE=[];
% alphaFastNewHipmean=[];
% alphaFastNewWOHip=[];
alphaSlowNew=[];
alphaSlowOldE=[];
% alphaSlowNewHipmean=[];
% alphaSlowNewWOHip=[];
% alphaSlowPlusbetaSlow=[];
% alphaSlowPlusbetaSlowmean=[];
% alphaSlowPlusbetaSlowWOHip=[];
% alphaFastPlusbetaFast=[];
% alphaFastPlusbetaFastmean=[];
% alphaFastPlusbetaFastWOHip=[];
% alphaSlowPlusXSlow=[];
% alphaSlowPlusXSlowmean=[];
% alphaSlowPlusXSlowWOHip=[];
% alphaFastPlusXFast=[];
% alphaFastPlusXFastmean=[];
% alphaFastPlusXFastWOHip=[];
% VcenterofoscillationNew=[];
% VcenterofoscillationmeanNew=[];
% VcenterofoscillationWOHipNew=[];
% VStepasymNew=[];
% VStepasymNewmean=[];
% VStepasymNewWOHip=[];


for i=1:1:length(condition)

alphaSlow=adaptData.getParamInCond('alphaTemp',condition(i));
% alphaSlowHipmean=adaptData.getParamInCond('alphaTempHipmean',condition(i));
% alphaSlowWOHip=adaptData.getParamInCond('alphaSlowWOHip',condition(i));


alphaFast=adaptData.getParamInCond('alphaFast',condition(i));
% alphaFastHipmean=adaptData.getParamInCond('alphaFastHipmean',condition(i));
% alphaFastWOHip=adaptData.getParamInCond('alphaFastWOHip',condition(i));


betaSlow=adaptData.getParamInCond('betaSlow',condition(i));
% betaSlowHipmean=adaptData.getParamInCond('betaSlowHipmean',condition(i));
% betaSlowWOHip=adaptData.getParamInCond('betaSlowWOHip',condition(i));

betaFast=adaptData.getParamInCond('betaFast',condition(i));
% betaFastHipmean=adaptData.getParamInCond('betaFastHipmean',condition(i));
% betaFastWOHip=adaptData.getParamInCond('betaFastWOHip',condition(i));

xSlow=adaptData.getParamInCond('xSlow',condition(i));
% xSlowHipmean=adaptData.getParamInCond('xSlowHipmean',condition(i));
% xSlowWOHip=adaptData.getParamInCond('xSlowWOHip',condition(i));

xFast=adaptData.getParamInCond('xFast',condition(i));
% xFastHipmean=adaptData.getParamInCond('xFastHipmean',condition(i));
% xFastWOHip=adaptData.getParamInCond('xFastWOHip',condition(i));

alphaFastN=alphaFast;
alphaFastOE=alphaFast;
% alphaFastNHipmean=alphaFastHipmean;
% alphaFastNWOHip=alphaFastWOHip;

alphaSlowN=alphaSlow;
alphaSlowOE=alphaSlow;
% % alphaSlowNHipmean=alphaSlowHipmean;
% % alphaSlowNWOHip=alphaSlowWOHip;

if strcmp(condition{i},'TM base')
RSlow=adaptData.getParamInCond('RSlowPos',condition(i));
RFast=adaptData.getParamInCond('RFastPos',condition(i));
RSlow=nanmean(RSlow);
RFast=nanmean(RFast);
end

if strcmp(condition{i},'Gradual Adaptation') || strcmp(condition{i},'Re-adaptation') || strcmp(condition{i},'adaptation') || strcmp(condition{i},'re-adaptation')
    
    for n=2:1:length(alphaSlow)
alphaFastN(n+1)=alphaSlowN(n)+betaSlow(n)-betaFast(n);
%         alphaFastNHipmean(n+1)=alphaSlowNHipmean(n)+betaSlowHipmean(n)-betaFastHipmean(n);
%         alphaFastNWOHip(n+1)=alphaSlowNWOHip(n)+betaSlowWOHip(n)-betaFastWOHip(n);
alphaFastOE(n)=alphaFastOE(n-1)+(1/(1+RFast))*(abs(betaFast(n-1))-(RFast*alphaFastOE(n-1)));
alphaSlowOE(n)=alphaSlowOE(n-1)+(1/(1+RSlow))*(abs(betaSlow(n-1))-(RSlow*alphaSlowOE(n-1)));

alphaSlowN(n)=alphaFastN(n-1)-xSlow(n-1)+xFast(n-1);
%         alphaSlowNHipmean(n)=alphaFastNHipmean(n-1)-xSlowHipmean(n-1)+xFastHipmean(n-1);
%         alphaSlowNWOHip(n)=alphaFastNWOHip(n-1)-xSlowWOHip(n-1)+xFastWOHip(n-1);       
    end
    
    alphaFastN=alphaFastN(1:end-1);
%     alphaFastNHipmean=alphaFastNHipmean(1:end-1);
%     alphaFastNWOHip=alphaFastNWOHip(1:end-1); 
end


% t=0;
% alphaSlowPlusbetaSlow2=[];
% alphaSlowPlusbetaSlowmean2=[];
% alphaSlowPlusbetaSlowWOHip2=[];
% alphaFastPlusbetaFast2=[];
% alphaFastPlusbetaFastmean2=[];
% alphaFastPlusbetaFastWOHip2=[];
% alphaSlowPlusXSlow2=[];
% alphaSlowPlusXSlowmean2=[];
% alphaSlowPlusXSlowWOHip2=[];
% alphaFastPlusXFast2=[];
% alphaFastPlusXFastmean2=[];
% alphaFastPlusXFastWOHip2=[];
% VcenterofoscillationN=[];
% VcenterofoscillationmeanN=[];
% VcenterofoscillationWOHipN=[];
% VStepasymNewN=[];
% VStepasymNewmeanN=[];
% VStepasymNewWOHipN=[];

%  for t=1:length(alphaSlowN)
           
%     alphaSlowPlusbetaSlow2(t,1)= (alphaSlowN(t)+betaSlow(t))/abs(alphaSlowN(t)-betaSlow(t));
%     alphaSlowPlusbetaSlowmean2(t,1)= (alphaSlowNHipmean(t)+betaSlowHipmean(t))/abs(alphaSlowNHipmean(t)-betaSlowHipmean(t));
%     alphaSlowPlusbetaSlowWOHip2(t,1)=(((betaSlowWOHip(t)-alphaSlowNWOHip(t))/2)+alphaSlowNWOHip(t))/abs(betaFastWOHip(t));
%     
%     alphaFastPlusbetaFast2(t,1)=(alphaFastN(t)+betaFast(t))/abs(alphaFastN(t)-betaFast(t));
%     alphaFastPlusbetaFastmean2(t,1)=(alphaFastNHipmean(t)+betaFastHipmean(t))/abs(alphaFastNHipmean(t)-betaFastHipmean(t));
%     alphaFastPlusbetaFastWOHip2(t,1)=(((betaFastWOHip(t)-alphaFastNWOHip(t))/2)+alphaFastNWOHip(t))/abs(betaFastWOHip(t));
%     
%     alphaSlowPlusXSlow2(t,1)=(alphaSlowN(t)+xSlow(t))/abs(alphaSlowN(t)-xSlow(t));
%     alphaSlowPlusXSlowmean2(t,1)=(alphaSlowNHipmean(t)+xSlowHipmean(t))/abs(alphaSlowNHipmean(t)-betaSlowHipmean(t));
%     alphaSlowPlusXSlowWOHip2(t,1)=(((xSlowWOHip(t)-alphaSlowNWOHip(t))/2)+alphaSlowNWOHip(t))/abs(betaFastWOHip(t));
%     
%     alphaFastPlusXFast2(t,1)=(alphaFastN(t)+xFast(t))/abs(alphaFastN(t)-xFast(t));
%     alphaFastPlusXFastmean2(t,1)=(alphaFastNHipmean(t)+xFastHipmean(t))/abs(alphaFastNHipmean(t)-betaFastHipmean(t));
%     alphaFastPlusXFastWOHip2(t,1)=(((xFastWOHip(t)-alphaFastNWOHip(t))/2)+alphaFastNWOHip(t))/abs(betaFastWOHip(t));
%     
%     VcenterofoscillationN(t,1)=alphaSlowPlusbetaSlow2(t,1)-alphaFastPlusbetaFast2(t,1);
%     VcenterofoscillationmeanN(t,1)=alphaSlowPlusbetaSlowmean2(t,1)-alphaFastPlusbetaFastmean2(t,1);
%     VcenterofoscillationWOHipN(t,1)=alphaSlowPlusbetaSlowWOHip2(t,1)-alphaFastPlusbetaFastWOHip2(t,1);
%     VStepasymNewN(t,1)=alphaSlowPlusXSlow2(t,1)-alphaFastPlusXFast2(t,1);
%     VStepasymNewmeanN(t,1)=alphaSlowPlusXSlowmean2(t,1)-alphaFastPlusXFastmean2(t,1);
%     VStepasymNewWOHipN(t,1)=alphaSlowPlusXSlowWOHip2(t,1)-alphaFastPlusXFastWOHip2(t,1);
%  end

alphaFastNew=[alphaFastNew; alphaFastN];
alphaFastOldE=[alphaFastOldE;alphaFastOE];

% alphaFastNewHipmean=[alphaFastNewHipmean; alphaFastNHipmean];
% alphaFastNewWOHip=[alphaFastNewWOHip; alphaFastNWOHip];

alphaSlowNew=[alphaSlowNew; alphaSlowN];
alphaSlowOldE=[alphaSlowOldE;alphaSlowOE];
% alphaSlowNewHipmean=[alphaSlowNewHipmean; alphaSlowNHipmean];
% alphaSlowNewWOHip=[alphaSlowNewWOHip; alphaSlowNWOHip];

% alphaSlowPlusbetaSlow=[alphaSlowPlusbetaSlow;alphaSlowPlusbetaSlow2];
% alphaSlowPlusbetaSlowmean=[alphaSlowPlusbetaSlowmean;alphaSlowPlusbetaSlowmean2];
% alphaSlowPlusbetaSlowWOHip=[alphaSlowPlusbetaSlowWOHip;alphaSlowPlusbetaSlowWOHip2];

% alphaFastPlusbetaFast=[alphaFastPlusbetaFast;alphaFastPlusbetaFast2];
% alphaFastPlusbetaFastmean=[alphaFastPlusbetaFastmean;alphaFastPlusbetaFastmean2];
% alphaFastPlusbetaFastWOHip=[alphaFastPlusbetaFastWOHip; alphaFastPlusbetaFastWOHip2];

% alphaSlowPlusXSlow=[alphaSlowPlusXSlow;alphaSlowPlusXSlow2];
% alphaSlowPlusXSlowmean=[alphaSlowPlusXSlowmean;alphaSlowPlusXSlowmean2];
% alphaSlowPlusXSlowWOHip=[alphaSlowPlusXSlowWOHip;alphaSlowPlusXSlowWOHip2];
% 
% alphaFastPlusXFast=[alphaFastPlusXFast;alphaFastPlusXFast2];
% alphaFastPlusXFastmean=[alphaFastPlusXFastmean;alphaFastPlusXFastmean2];
% alphaFastPlusXFastWOHip=[alphaFastPlusXFastWOHip;alphaFastPlusXFastWOHip2];

% VcenterofoscillationNew=[VcenterofoscillationNew;VcenterofoscillationN];
% VcenterofoscillationmeanNew=[VcenterofoscillationmeanNew;VcenterofoscillationmeanN];
% VcenterofoscillationWOHipNew=[VcenterofoscillationWOHipNew;VcenterofoscillationWOHipN];
% VStepasymNew=[VStepasymNew;VStepasymNewN];
% VStepasymNewmean=[VStepasymNewmean;VStepasymNewmeanN];
% VStepasymNewWOHip=[VStepasymNewWOHip;VStepasymNewWOHipN];


end


% this=paramData([adaptData.data.Data,alphaFastNew,alphaFastNewHipmean,alphaFastNewWOHip,alphaSlowNew,alphaSlowNewHipmean,alphaSlowNewWOHip,...
%     alphaSlowPlusbetaSlow,alphaSlowPlusbetaSlowmean,alphaSlowPlusbetaSlowWOHip,alphaFastPlusbetaFast,alphaFastPlusbetaFastmean,...
%     alphaFastPlusbetaFastWOHip,alphaSlowPlusXSlow,alphaSlowPlusXSlowmean,alphaSlowPlusXSlowWOHip,alphaFastPlusXFast,alphaFastPlusXFastmean,alphaFastPlusXFastWOHip,...
%     VcenterofoscillationNew,VcenterofoscillationmeanNew,VcenterofoscillationWOHipNew,VStepasymNew,VStepasymNewmean,VStepasymNewWOHip],...
%     [adaptData.data.labels 'alphaFastNew' 'alphaFastNewHipmean' 'alphaFastNewWOHip' 'alphaSlowNew' 'alphaSlowNewHipmean'...
%     'alphaSlowNewWOHip' 'alphaSlowPlusbetaSlowNew' 'alphaSlowPlusbetaSlowmeanNew' 'alphaSlowPlusbetaSlowWOHipNew' 'alphaFastPlusbetaFastNew' 'alphaFastPlusbetaFastmeanNew'...
%     'alphaFastPlusbetaFastWOHipNew' 'alphaSlowPlusXSlowNew' 'alphaSlowPlusXSlowmeanNew' 'alphaSlowPlusXSlowWOHipNew' 'alphaFastPlusXFastNew' 'alphaFastPlusXFastmeanNew' 'alphaFastPlusXFastWOHipNew'...
%  'VcenterofoscillationNew' 'VcenterofoscillationmeanNew' 'VcenterofoscillationWOHipNew' 'VStepasymNew' 'VStepasymNewmean' 'VStepasymNewWOHip'],adaptData.data.indsInTrial,adaptData.data.trialTypes);
% adaptData=adaptationData(rawExpData.metaData,rawExpData.subData,this);
 this=paramData([adaptData.data.Data,alphaFastNew,alphaSlowNew,alphaFastOldE,alphaSlowOldE],[adaptData.data.labels 'alphaFastNew' 'alphaSlowNew' 'alphaFastOLD' 'alphaSlowOLD'],adaptData.data.indsInTrial,adaptData.data.trialTypes);
 adaptData=adaptationData(rawExpData.metaData,rawExpData.subData,this);
 
 saveloc=[];
%  save([saveloc subject 'simulation.mat'],'adaptData'); 
 cd('C:\Users\dum5\Desktop\dulce\Exp0002\matrix')
 save([saveloc subject 'simulation.mat'],'adaptData'); 