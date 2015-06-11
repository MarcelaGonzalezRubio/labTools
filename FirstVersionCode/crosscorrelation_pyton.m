load('NBF15.mat')

NexuslowFreq=[];
RHSlowFreq=[];

GRRz=expData.data{4}.GRFData.Data(:,9);
GRLz=expData.data{4}.GRFData.Data(:,10);
for i=1:10:length(GRRz)  
    NexuslowFreq1=GRRz(i);
    NexuslowFreq=[NexuslowFreq;NexuslowFreq1];
end

[header,data] = JSONtxt2cell('1423062582data.txt');

newData=nan((((data(end,1)-data(1,1)))+1),9);
newData2=nan((((data(end,1)-data(1,1)))+1),9);


data(:,1)=data(:,1)-data(1,1)+1;

gap=diff(data(:,1));
figure()
plot(gap,'b')

newData=interp1(data(:,1),data(:,1:end),[data(1,1):data(end,1)],'linear');

for i=1:length(data)
        newData2(data(i,1),1:end)=data(i,1:end);
end

GRRz=NexuslowFreq; 
FzR=newData(:,2);

%Trial 1
[acor, lag]=xcorr(GRRz,FzR);
[~,I]=max((acor));
timeDiff=lag(I);
 B=[];
newData=newData(abs(timeDiff)+1:end,1:end);
newData2=newData2(abs(timeDiff)+1:end,1:end);
 
 
figure()
plot(GRRz,'b')
hold on 
plot(newData(:,2), 'r')
plot(newData2(:,2), 'g')

newData=newData(abs(timeDiff)+1:end,1:end);
newData2=newData2(abs(timeDiff)+1:end,1:end);

%Detection of HS for R leg with nexus data 
% for z=1:length(NexuslowFreq)-10
%     if NexuslowFreq(z)>-30 && NexuslowFreq(z+1)<-30
%         RHSlowFreq(z)=1;
%     else
%         RHSlowFreq(z)=0;
%     end
% end 
% % 
z=0;
%Detection of HS 
%  for z=1:length(newData)
%      % R leg
%     if newData(z,2)>-30 && newData(z+1,2)<=-30
%         newRHSlowFreq(z,1)=1;
%     else
%         newRHSlowFreq(z,1)=0;
%     end
%     
%     %L leg
%     if newData(z,3)>-30 && newData(z+1,3)<=-30
%         newLHSlowFreq(z,1)=1;
%     else
%         newLHSlowFreq(z,1)=0;
%     end
%     
%  end
function [LHSnexus,RHSnexus]= getEventsFromForces(,FzR,fsample)
function [LHSpyton,RHSpyton]= getEventsFromForces(newData(:,3),newData(:,2),100)
% Index R leg 
 index=find(newRHSlowFreq);
 
 for n=1:length(index)-1
     if index(n+1) -index(n)<100
      newRHSlowFreq(index(n))=0;
     end
 end
 
 n=0;
 % Index L leg 
  index2=find(newLHSlowFreq);
  for n=1:length(index2)-1
     if index2(n+1) -index2(n)<100
      newLHSlowFreq(index2(n))=0;
     end
 end
 % GOOD STRIDES
 for x=1:length(newData2)
    if newRHSlowFreq(x)==1 && newData2(x+10,4)==1
        GoodRHS(x,1)=1;
    else 
        GoodRHS(x,1)=0;
    end
    
     if newLHSlowFreq(x)==1 && newData2(x+10,5)==1
        GoodLHS(x,1)=1;
    else 
        GoodLHS(x,1)=0;
    end 
 end
 
 TotalGoodRHS=sum(GoodRHS);
 TotalGoodLHS=sum(GoodLHS);
 TotalRHS=sum(newRHSlowFreq);
 TotalLHS=sum(newLHSlowFreq);
 
 PGoodRHS=TotalGoodRHS/TotalRHS;
 PGoodLHS=TotalGoodLHS/TotalLHS;
 
 figure()
 for i=1:1:1
hold on
bar((1:1)+(.5+.5.*i),PGoodRHS(i,1),0.2,'FaceColor',[.8,.8,.8])
bar((1:1)+(.7+.5*i),PGoodLHS(i,1),0.2,'FaceColor',[.0,.36,.6])
 end

condition={'Gradual Adaptation','Re adaptation'};
xTickPos=2:.5:2*length(condition);
set(gca,'XTick',xTickPos,'XTickLabel',condition)
legend( 'Fast Leg','Slow Leg')
%title(['Good Steps' '(',Subject ')'])
%%
figure()
plot(newRHSlowFreq*100) 
hold on
plot(GoodRHS*100,'mo','MarkerSize',5)
plot(newLHSlowFreq*100,'r') 
plot(GoodLHS*100,'gx','MarkerSize',5)