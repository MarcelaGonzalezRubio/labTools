% function a=crosscorrelation_pyton_rev3(subject,dynamic)
% for w=1:2
subject=('PDT04');
dynamic=1;
load([subject 'params.mat'])
load([subject '.mat'])
load([subject 'RAW.mat'])
w=1;

load(['Pyton' num2str(w) '.mat'])

j(1)=adaptData.metaData.trialsInCondition{3};
% j(2)=adaptData.metaData.trialsInCondition{5};


%Force plate data from Nexus
GRRz=expData.data{j(w)}.GRFData.Data(:,9);
GRLz=expData.data{j(w)}.GRFData.Data(:,3);

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


%Determination of crosscorrelation between Nexus at 100Hz and Interpolate
%Pyton data
[acor, lag]=xcorr(NexusRlowFreq,newData(:,2));
[~,I]=max((acor));
timeDiff=lag(I);
newData=newData(abs(timeDiff)+1:end,1:end);
newData2=newData2(abs(timeDiff)+1:end,1:end);
% 
NexusLlowFreq=NexusLlowFreq(1:length(newData));
NexusRlowFreq=NexusRlowFreq(1:length(newData)); 

% figure()
% plot(NexusRlowFreq,'b')
% hold on 
% plot(newData(:,2), 'r')
% plot(newData2(:,2), 'g')
% legend('Nexus','Interpolate Pyton','Pyton')
% title('Sync of R Force plate data')

% interpolateData=labTimeSeries(newData,0,0.01,{'FrameNumber','Rfz','Lfz','RHS','LHS','RGOOD','LGOOD','Ralpha','Lalpha','R target','L target','RHIP','LHIP','RANK','LANK'});
% NaNdata=labTimeSeries(newData2,0,0.01,{'FrameNumber','Rfz','Lfz','RHS','LHS','RGOOD','LGOOD','Ralpha','Lalpha','R target','L target','RHIP','LHIP','RANK','LANK'});
%%
%Finding HS from Nexus at 100HZ and Interpolated Pyton data, interpolate
%data is used to make sure that we dont take in consideration extras HS.
[LHSnexus,RHSnexus]= getEventsFromForces(NexusLlowFreq,NexusRlowFreq,100);
[LHSpyton,RHSpyton]= getEventsFromForces(newData(:,3),newData(:,2),100);

%% 
%localication of HS
locRHSpyton=find(RHSpyton==1);
locRHSnexus=find(RHSnexus==1);
locLHSpyton=find(LHSpyton==1);
locLHSnexus=find(LHSnexus==1);
% locRindex=locRHSpyton;
% locLindex=locLHSpyton;
locRindex=find(newData(:,4)==1);
locLindex=find(newData(:,5)==1);

%Delete extras HS deteted by Python 
while length(locRHSpyton)~=length(locRindex)
    diffLengthR=length(locRindex)-length(locRHSpyton);
   FrameDiffR=locRindex(1:end-diffLengthR)-locRHSpyton;
   IsBadR=find(FrameDiffR<=-10);
   if isempty(IsBadR)
       stop
       display('Something is WRONG!')
   else
       locRindex(IsBadR(1))=[];
   end
end

while length(locLHSpyton)~=length(locLindex)
    diffLength=length(locLindex)-length(locLHSpyton);
   FrameDiff=locLindex(1:end-diffLength)-locLHSpyton;
   IsBad=find(FrameDiff<=-10);
   if isempty(IsBad)
       stop
       display('Something is WRONG!')
   else
       locLindex(IsBad(1))=[];
   end
end

%%
%Good strides
%GoodEvents1=expData.data{4}.adaptParams.Data(:,1);
GoodEvents=expData.data{j(w)}.adaptParams.Data(:,1); %Good events detected
GoodRHS=newData(locRindex,6);
GoodLHS=newData(locLindex,7); 
GoodRHS=GoodRHS(GoodEvents==1);
GoodLHS=GoodLHS(GoodEvents==1);
%hola=nan(length(GoodEvents1),1);
% hola=[];
% GoodRHS2=[hola;GoodRHS];
% GoodLHS2=[hola;GoodLHS];


%%
%find alpha value on time
alphaR_time=nan(length(newData),1);
alphaL_time=nan(length(newData),1);
alphaR_time(locRindex,1)=newData(locRindex,8)*1000;
alphaL_time(locLindex,1)=newData(locLindex,9)*1000;
%alpha values at HS 
alphaRPyton=newData(locRindex,8)*1000;
alphaLPyton=newData(locLindex,9)*1000;
alphaRPytonGood=alphaRPyton(GoodEvents==1);
alphaLPytonGood=alphaLPyton(GoodEvents==1);
% 
if dynamic ==1 
Rtarget=newData(locRindex,10)*1000;
Ltarget=newData(locLindex,11)*1000;
RtargetGood=Rtarget(GoodEvents==1);
LtargetGood=Ltarget(GoodEvents==1);
elseif dynamic== 0 %static target
Rscale=newData(locRindex,10);
Lscale=newData(locLindex,11);
RscaleGood=Rscale(GoodEvents==1);
LscaleGood=Lscale(GoodEvents==1);
RtargetGood=(0.5./RscaleGood)*1000;
LtargetGood=(0.5./LscaleGood)*1000;    
end

%%
%RAW data 
 
% RHIP=expData.data{j(w)}.markerData.Data(:,14);
% RANK=expData.data{j(w)}.markerData.Data(:,26);
% LHIP=expData.data{j(w)}.markerData.Data(:,29);
% LANK=expData.data{j(w)}.markerData.Data(:,47);
% Rleg=RHIP-RANK;
% Lleg=LHIP-LANK;
% 
% alphaRNexus=RHIP(locRHSnexus)-RANK(locRHSnexus);
% alphaLNexus=LHIP(locLHSnexus)-LANK(locLHSnexus);
% alphaRNexus_time=nan(length(Rleg),1);
% alphaLNexus_time=nan(length(Rleg),1);
% alphaRNexus_time(locRHSnexus,1)=RHIP(locRHSnexus)-RANK(locRHSnexus);
% alphaLNexus_time(locLHSnexus,1)=LHIP(locLHSnexus)-LANK(locLHSnexus);
% 
% alphaRnexusGood=alphaRNexus(GoodEvents==1);
% alphaLnexusGood=alphaLNexus(GoodEvents==1);
% this=paramData([adaptData.data.Data,GoodRHS2,GoodLHS2],[adaptData.data.labels 'GoodStrideR' 'GoodStrideL'],adaptData.data.indsInTrial,adaptData.data.trialTypes);
% adaptData=adaptationData(rawExpData.metaData,rawExpData.subData,this);
% saveloc=[];
% save([saveloc subject 'params.mat'],'adaptData'); 

%%       
%  TotalGoodRHS=sum(GoodRHS==1);
%  TotalGoodLHS=sum(GoodLHS==1);
%  TotalRHS=sum(RHSpyton);
%  TotalLHS=sum(LHSpyton);
%  
%  PGoodRHS=TotalGoodRHS/TotalRHS;
%  PGoodLHS=TotalGoodLHS/TotalLHS;
%  
%  figure()
%  for i=1:1:1
% hold on
% bar((1:1)+(.5+.5.*i),PGoodRHS(i,1),0.2,'FaceColor',[.8,.8,.8])
% bar((1:1)+(.7+.5*i),PGoodLHS(i,1),0.2,'FaceColor',[.0,.36,.6])
%  end
% 
% condition={'Gradual Adaptation','Re adaptation'};
% xTickPos=2:.5:2*length(condition);
% set(gca,'XTick',xTickPos,'XTickLabel',condition)
% legend( 'Fast Leg','Slow Leg')
% title(['Good Steps' '(',subject ')'])
%%
% figure()
% plot(RHSpyton*100) 
% hold on
% plot(GoodRHS*100,'mo','MarkerSize',5)
% plot(LHSpyton*100,'r') 
% plot(GoodLHS*100,'gx','MarkerSize',5)
% % end
%%
%Comporaring RAW data with pyton recording

% figure()
% plot(RHIP)
% hold on
% plot(newData(:,12)*1000,'r')
% plot(newData2(:,12)*1000,'g')
% title('RHip marker position')
% legend('Nexus','Interpolate Pyton','Pyton')
% % 
% figure()
% plot(LHIP)
% hold on
% plot(newData(:,13)*1000,'r')
% plot(newData2(:,13)*1000,'g')
% title('LHip marker position')
% legend('Nexus','Interpolate Pyton','Pyton')
% 
% figure()
% plot(RANK)
% hold on
% plot(newData(:,14)*1000,'r')
% plot(newData2(:,14)*1000,'g')
% title('RANK marker position')
% legend('Nexus','Interpolate Pyton','Pyton')
% 
% figure()
% plot(LANK)
% hold on
% plot(newData(:,15)*1000,'r')
% plot(newData2(:,15)*1000,'g')
% title('LANK marker position')
% legend('Nexus','Interpolate Pyton','Pyton')

 %%
% figure()
% hold on
% plot(Rleg,'b')
% plot(newData(:,8)*1000,'r')
% plot(alphaRNexus_time,'og','MarkerSize',3)
% plot(alphaR_time,'ok','MarkerSize',3)
% legend('Right leg position Nexus','Right leg position Pyton','RHS nexus','RHS Pyton')
% title(['RHIP-RANK in time with detection of HS' '(',subject ')'])
% 
% figure()
% hold on
% plot(Lleg,'b')
% plot(newData(:,9)*1000,'r')
% plot(alphaLNexus_time,'og','MarkerSize',3)
% plot(alphaL_time,'ok','MarkerSize',3)
% legend('Left leg position Nexus','Left leg position Pyton','LHS nexus','LHS Pyton')
% title(['LHIP-LANK in time with detection of HS'  '(',subject ')'])
%%
% alphaRpost=adaptData.getParamInCond('alphaFast','Gradual Adaptation');
% alphaLpost=adaptData.getParamInCond('alphaSlow','Gradual Adaptation');


%%
% if dynamic==0
% ystd=25*ones([length(GoodRHS),1]);
% 
% else
ystdR=25*ones([length(GoodRHS),1]);
ystdL=25*ones([length(GoodRHS),1]);
% end 
figure()
hold on 

errorbar(RtargetGood,ystdR,'rx');

for i=1:length(GoodRHS)
        %if GoodRHS(i)==0
        if abs(alphaRPytonGood(i)-RtargetGood(i))<25
         plot(i,alphaRPytonGood(i),'g.','MarkerSize',15)
          GoodR(i,1)=1;
    %elseif GoodRHS(i)==1
        elseif abs(alphaRPytonGood(i)-RtargetGood(i))>=25
        plot(i,alphaRPytonGood(i),'k.','MarkerSize',15)
         GoodR(i,1)=0;
        end
    
end
title(['Alpha R leg' '(',subject ')' '   \alpha_(_n_)-\alpha^*<25' ])
axis tight
figure()
hold on 
errorbar(LtargetGood,ystdL,'rx');

for i=1:length(GoodLHS)
    %if GoodLHS(i)==0
    if abs(alphaLPytonGood(i)-LtargetGood(i))<25
        plot(i,alphaLPytonGood(i),'g.','MarkerSize',15)
         GoodL(i,1)=1;
    %elseif GoodLHS(i)==1
     elseif abs(alphaLPytonGood(i)-LtargetGood(i))>=25
        plot(i,alphaLPytonGood(i),'k.','MarkerSize',15)
         GoodL(i,1)=0;
    end
  
end
title(['Alpha L leg' '(' ,subject ')' '  \alpha_(_n_)-\alpha^*<25' ])
axis tight
%legend('Bad Steps','God Steps')
%%
% figure()
% plot(NexusRlowFreq,'b')
% hold on 
% plot(newData(:,2), 'r')
% plot(RHSpyton*100, 'm')
% plot(RHSnexus*100,'--k')
% plot(newData(:,4)*100,'g')
% title('Ground reaction with HS detection ')
% legend('Nexus GRz','Pyton GRz','RHS pyton using Groud Reactio','RHS nexus','RHSpyton interpolate')
% end
% end
%%
% creating vector for good or bad steps 
condition= adaptData.metaData.conditionName;
% This is new ~~~~~~~~~~~~~~~~~
     %Needed in case subject did not perform one of the conditions
        %in the condition list
        
condition=condition(find(~cellfun(@isempty,adaptData.metaData.trialsInCondition)));
StepsR=[];
StepsL=[];
for i=1:length(condition)
    
    numberSteps=adaptData.getParamInCond('Good',condition(i));
    
    StepsR2=NaN(length(numberSteps),1);
    StepsL2=NaN(length(numberSteps),1);
    
    if strcmp(condition{i},'Gradual adaptation') 
    StepsR2=GoodRHS;
    StepsL2=GoodLHS;
    end
    
  StepsR=[StepsR;StepsR2];
  StepsL=[StepsL;StepsL2];
end
    
    
    