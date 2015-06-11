function [results]=meanTrials(subject,removeBiasFlag)

load([subject 'params.mat'])

adaptData.removeBias;
condition= ans.metaData.conditionName;
labels=ans.data.labels;

results.avg=[];
results.std=[];
results.rms=[];

for i=1:1:length(condition)

    StepLength=[];
 
%Fast leg   
StepLength=ans.getParamInCond('stepLengthAsym',condition(i));

results.avg.fast(1,i)=nanmean(StepLength);
results.std.fast(1,i)=nanstd(StepLength);
results.rms.fast(1,i)=rms(StepLength);
results.stderror.fast(1,i)=(results.std.fast(1,i))/(sqrt(length(StepLength)));

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


%BAR PLOT FOR AVG AND STD ERROR
figure()
for i=1:1:length(condition)

hold on
bar((1:1)+(.5+.5.*i),Newresults.avg.fast(1,i),0.2,'FaceColor',[.8,.8,.8])
%bar((1:1)+(.7+.5*i),Newresults.avg.slow(1,i),0.2,'FaceColor',[.0,.36,.6])
errorbar(((1:1)+(.5+.5*i)),Newresults.avg.fast(1,i),Newresults.stderror.fast(1,i),'.r')
%errorbar(((1:1)+(.7+.5*i)),Newresults.avg.slow(1,i),Newresults.stderror.slow(1,i),'.r')
end

title(['Avg per Trial' ' (',subject ')'])
xTickPos=2:.5:2*length(condition);
set(gca,'XTick',xTickPos,'XTickLabel',condition)
%legend( 'Fast Leg','Slow Leg')
hold off

end
