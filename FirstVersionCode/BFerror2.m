function [results]=BFerror2(subject,condition,removeBiasFlag) 

load([subject 'params.mat'])

EFastT=[];
ESlowT=[];
results.avg=[];
results.std=[];
results.rms=[];

for i=1:1:length(condition)

EFast=[];
ESlow=[];
betaFast=[];
alphaFast=[];
RatioFast=[];
betaSlow=[];
alphaSlow=[];
RatioSlow=[];

%Fast leg   
betaFast=adaptData.getParamInCond('betaFast',condition(i));
alphaFast=adaptData.getParamInCond('alphaFast',condition(i));
RatioFast=adaptData.getParamInCond('RFastPos',condition(i));
 
if length(RatioFast)>=290
     RatioFast=RatioFast([1:290],1);
else
    RatioFast=RatioFast;
end 


%Slow Leg
betaSlow=adaptData.getParamInCond('betaSlow',condition(i));
alphaSlow=adaptData.getParamInCond('alphaTemp',condition(i));
RatioSlow=adaptData.getParamInCond('RSlowPos',condition(i));

if length(RatioSlow)>=290
RatioSlow=RatioSlow([1:290],1);
else
    RatioSlow=RatioSlow;
end

avgFast=median(RatioFast);
avgSlow=median(RatioSlow);

%Error calculation e=|beta|-|R*alpha|
for n=1:1:length(betaFast)
    
    EFast(n,1)=abs(betaFast(n,1))-abs(avgFast*alphaFast(n,1));
    ESlow(n,1)=abs(betaSlow(n,1))-abs(avgSlow*alphaSlow(n,1));
    
end

results.avg.fast(1,i)=nanmean(EFast);
results.std.fast(1,i)=nanstd(EFast);
results.rms.fast(1,i)=rms(EFast);
results.stderror.fast(1,i)=(results.std.fast(1,i))/(sqrt(length(EFast)));

results.avg.slow(1,i)=nanmean(ESlow);
results.std.slow(1,i)=nanstd(ESlow);
results.rms.slow(1,i)=rms(ESlow);
results.stderror.slow(1,i)=(results.std.slow(1,i))/(sqrt(length(ESlow)));

EFastT=[EFastT;EFast];  
ESlowT=[ESlowT;ESlow]; 



end

if removeBiasFlag==1
  %[a,b,c]=intersect('TM base', condition);
  %for z=1:1:length(condition)
      %Newresults.avg.fast(1,z)=results.avg.fast(1,z)-results.avg.fast(1,c);
      %Newresults.avg.slow(1,z)=results.avg.slow(1,z)-results.avg.slow(1,c);
  %end
   Newresults=results;
else
    Newresults=results;
end

% 
% figure()
% subplot(2,1,1),plot(EFastT,'.','MarkerSize',15)
% title(subject)
% ylabel('Error Fast')
% axis([0 2000 -200 200])
% subplot(2,1,2),plot(ESlowT,'.','MarkerSize',15)
% title(subject)
% ylabel('Error Slow')
% axis([0 2000 -200 200])

%BAR PLOT FOR AVG AND STD ERROR
figure()
for i=1:1:length(condition)

hold on
bar((1:1)+(.5+.5.*i),Newresults.avg.fast(1,i),0.2,'FaceColor',[.8,.8,.8])
bar((1:1)+(.7+.5*i),Newresults.avg.slow(1,i),0.2,'FaceColor',[.0,.36,.6])
errorbar(((1:1)+(.5+.5*i)),Newresults.avg.fast(1,i),Newresults.stderror.fast(1,i),'.r')
errorbar(((1:1)+(.7+.5*i)),Newresults.avg.slow(1,i),Newresults.stderror.slow(1,i),'.r')
end

title(['Avg Error' ' (',subject ')'])
xTickPos=2.1:.5:2*length(condition);
set(gca,'XTick',xTickPos,'XTickLabel',condition)
legend( 'Fast Leg','Slow Leg')
hold off

end
