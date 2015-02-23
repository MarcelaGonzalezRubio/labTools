%function crosscorrelation_pyton_rev2(subject)
subject=('PDT02');
load([subject 'params.mat'])
load([subject '.mat'])
load('Pyton.mat');

NexusRlowFreq=[];
NexusLlowFreq=[];

j=adaptData.metaData.trialsInCondition{3};

%Force plate data from Nexus
GRRz=expData.data{j}.GRFData.Data(:,9);
GRLz=expData.data{j}.GRFData.Data(:,3);

%Converting force plate data from 1000Hz to 100Hz
NexusRlowFreq=resample(GRRz,1,10);
NexusLlowFreq=resample(GRLz,1,10);


%Creating NaN matrix with the lenght of the data 
newData=nan((((outmat(end,1)-outmat(1,1)))+1),15);
newData2=nan((((outmat(end,1)-outmat(1,1)))+1),15);

%Making frames from Pyton start at 1 
outmat(:,1)=outmat(:,1)-outmat(1,1)+1;

data=unique(outmat(:,1));

%finding unique colums
for i=1:length(data)
[datarows(i),~]=find(outmat(:,1)==data(i),1,'first');
end

outmat=outmat(datarows,:);
%Calculating the gap's length for pyton data
% gap=diff(outmat(:,1));
% figure()
% plot(gap,'b')

%Creating a linear interpolate matrix from Pyton data 
newData=interp1(outmat(:,1),outmat(:,1:end),[outmat(1,1):outmat(end,1)]);

%Creating a Matrix with NaN in gaps from Pyton 
for i=1:length(outmat)
        newData2(outmat(i,1),1:end)=outmat(i,:);
end


GRRz=NexusRlowFreq; 
FzR=newData(:,2);

%Determination of crosscorrelation between Nexus at 100Hz and Interpolate
%Pyton data
[acor, lag]=xcorr(GRRz,FzR);
[~,I]=max((acor));
timeDiff=lag(I);
newData=newData(abs(timeDiff)+1:end,1:end);
newData2=newData2(abs(timeDiff)+1:end,1:end);
 
%figure()
plot(GRRz,'b')
hold on 
plot(newData(:,2), 'r')
plot(newData2(:,2), 'g')
legend('Nexus','Interpolate Pyton','Pyton')
title('Sync of R Force plate data')

%%
%Finding HS from Nexus at 100HZ and Interpolated Pyton data
[LHSnexus,RHSnexus]= getEventsFromForces(NexusLlowFreq,NexusRlowFreq,100);
[LHSpyton,RHSpyton]= getEventsFromForces(newData(:,3),newData(:,2),100);

%%
%localication of HS
locRHSpyton=find(RHSpyton==1);
locRHSnexus=find(RHSnexus==1);
locLHSpyton=find(LHSpyton==1);
locLHSnexus=find(LHSnexus==1);

%%
%Good strides 
GoodRHS=newData(locRHSpyton+10,6);
GoodLHS=newData(locLHSpyton+10,6);
%%
% for x=1:length(locRHSpyton)
%     if RHSpyton(locRHSpyton(x))==1 && newData(locRHSpyton(x)+10,6)==1
%         GoodRHS(x,1)=1;
%         
%     elseif  RHSpyton(locRHSpyton(x))==1 && newData(locRHSpyton(x)+10,6)==0
%         GoodRHS(x,1)=0;
%     
%     end
%       
% end


% for x=1:length(RHSpyton)
%     if RHSpyton((x))==1 && newData((x)+10,6)==1
%         GoodRHS(x,1)=1;
%         
%     elseif  RHSpyton((x))==1 && newData((x)+10,6)==0
%         GoodRHS(x,1)=0;
%     
%     end
%       
% end

%  for x=1:length(locLHSpyton)
%  if LHSpyton(locLHSpyton(x))==1 && newData(locLHSpyton(x)+10,7)==1
%         GoodLHS(x,1)=1;
%         
%  elseif  LHSpyton(locLHSpyton(x))==1 && newData(locLHSpyton(x)+10,7)==0
%         GoodLHS(x,1)=0;
%          
%     end 
%  end
 
 %%
 alphaRPyton(:,1)=newData(locRHSpyton,8)*1000;
 alphaLPyton(:,1)=newData(locLHSpyton,9)*1000;
 
 %%
  for x=1:length(newData)
      if RHSpyton(x)==1 
          GOOD(x,1)=1;
        alphaTEST(x,1)=newData(x,8);
         
      end 
      if LHSpyton(x)==1 
          GOODLL(x,1)=1;
        alphaTESTL(x,1)=newData(x,9);
      end
  end
  
RHIP=expData.data{j}.markerData.Data(:,14);
RANK=expData.data{j}.markerData.Data(:,26);
LHIP=expData.data{j}.markerData.Data(:,29);
LANK=expData.data{j}.markerData.Data(:,47);
Rleg=RHIP-RANK;

alphaRNexus(:,1)=RHIP(locRHSpyton)-RANK(locRHSpyton);
alphaLNexus(:,1)=LHIP(locLHSpyton)-LANK(locLHSpyton);
       for x=1:length(RHIP)
           if RHSnexus(x)==1
               alphaRNExus(x,1)=RHIP(x)-RANK(x);
           else
                alphaRNExus(x,1)=0;
           end
           
       end
       
       
 %
 TotalGoodRHS=sum(GoodRHS);
 TotalGoodLHS=sum(GoodLHS);
 TotalRHS=sum(RHSpyton);
 TotalLHS=sum(LHSpyton);
 
 PGoodRHS=TotalGoodRHS/TotalRHS;
 PGoodLHS=TotalGoodLHS/TotalLHS;
 
 figure()
 for i=1:1:1
hold on
bar((1:1)+(.5+.5.*i),PGoodRHS(i,1),0.2,'FaceColor',[.8,.8,.8])
bar((1:1)+(.7+.5*i),PGoodLHS(i,1),0.2,'FaceColor',[.0,.36,.6])
 end
% 
% condition={'Gradual Adaptation','Re adaptation'};
% xTickPos=2:.5:2*length(condition);
% set(gca,'XTick',xTickPos,'XTickLabel',condition)
% legend( 'Fast Leg','Slow Leg')
% %title(['Good Steps' '(',Subject ')'])
%%
figure()
plot(RHSpyton*100) 
hold on
plot(GoodRHS*100,'mo','MarkerSize',5)
plot(LHSpyton*100,'r') 
plot(GoodLHS*100,'gx','MarkerSize',5)
% end
%%
figure()
plot(RHIP)
hold on
plot(newData(:,12)*1000,'r')
plot(newData2(:,12)*1000,'g')
title('RHip marker position')

figure()
plot(LHIP)
hold on
plot(newData(:,13)*1000,'r')
plot(newData2(:,13)*1000,'g')
title('LHip marker position')

figure()
plot(RANK)
hold on
plot(newData(:,14)*1000,'r')
plot(newData2(:,14)*1000,'g')
title('RANK marker position')


figure()
plot(LANK)
hold on
plot(newData(:,15)*1000,'r')
plot(newData2(:,15)*1000,'g')
title('LANK marker position')


