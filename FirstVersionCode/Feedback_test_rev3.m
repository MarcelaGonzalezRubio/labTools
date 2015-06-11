%Simulation of feedback 
function []=Feedback_test(subject)

 load([subject 'params.mat'])
 load([subject 'RAW.mat'])

condition= adaptData.metaData.conditionName;
% This is new ~~~~~~~~~~~~~~~~~
     %Needed in case subject did not perform one of the conditions
        %in the condition list
condition=condition(find(~cellfun(@isempty,adaptData.metaData.trialsInCondition)));


for i=1:1:length(condition)

alphaSlow=adaptData.getParamInCond('alphaTemp',condition(i));
alphaFast=adaptData.getParamInCond('alphaFast',condition(i));
TargetFast=[];
TargetSlow=[];
TargetFastOld=[];
TargetSlowOld=[];

if strcmp(condition{i},'TM base')
RSlow=adaptData.getParamInCond('RSlowPos',condition(i));
RFast=adaptData.getParamInCond('RFastPos',condition(i));
RSlow=nanmean(RSlow);
RFast=nanmean(RFast);
end



if strcmp(condition{i},'Gradual Adaptation') || strcmp(condition{i},'Re-adapation') || strcmp(condition{i},'adaptation')
n=[];
h=[];
TargetFast=[];
TargetSlow=[];
TargetFastOld=[];
TargetSlowOld=[];
GoodR=[];

    betaSlow=adaptData.getParamInCond('betaSlow',condition(i));
    betaFast=adaptData.getParamInCond('betaFast',condition(i));
    xSlow=adaptData.getParamInCond('xSlow',condition(i));
    xFast=adaptData.getParamInCond('xFast',condition(i));
    
    
    for n=[4:4:length(alphaSlow)]
        TargetFast(n+1:n+4)=alphaFast(n)+betaSlow(n)-betaFast(n);
        TargetSlow(n:n+3)=alphaSlow(n-1)-xSlow(n-1)+xFast(n-1);
        
        TargetFastOld(n:n+3)=alphaFast(n-1)+(1/(1+RFast))*(abs(betaFast(n-1))-RFast*(alphaFast(n-1)));
        TargetSlowOld(n:n+3)=alphaSlow(n-1)+(1/(1+RSlow))*(abs(betaSlow(n-1))-RSlow*(alphaSlow(n-1)));
    end
    
    TargetFast=TargetFast';
    TargetSlow=TargetSlow';
    TargetFastOld=TargetFastOld';
    TargetSlowOld=TargetSlowOld';
    
    for h=1:1:length(alphaFast)  
    
    if TargetFastOld(h)-25<= alphaFast(h) && alphaFast(h)<= TargetFastOld(h)+25
        GoodR(h,1)=1;
    else
        GoodR(h,1)=0;
    end
    
    if TargetSlowOld(h)-25<=alphaSlow(h) && alphaSlow(h)<=TargetSlowOld(h)+25
        GoodL(h,1)=1;
    else
        GoodL(h,1)=0;
    end
    end
    
end


alphaFastN=[alphaFast; TargetFast];
alphaFastO=[alphaFast; TargetFastOld];
alphaSlowN=[alphaSlow; TargetSlow];
alphaSlowO=[alphaSlow; TargetSlowOld];

end
