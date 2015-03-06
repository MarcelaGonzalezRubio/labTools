%SyncronizationPython 
 function a=syncronizationPython(subject,dynamic)

GRRz=[];
GRLz=[];

 
 load([subject 'params.mat'])
load([subject '.mat'])
load([subject 'RAW.mat'])
load(['Pyton' num2str(w) '.mat'])


j(w)=adaptData.metaData.trialsInCondition{3};

%Ground reaction data from Nexus
GRRz=expData.data{j(w)}.GRFData.Data(:,9);
GRLz=expData.data{j(w)}.GRFData.Data(:,3);

%resample GR from 1000Hz to 100Hz
NexusRlowFreq=resample(GRRz,1,10);
NexusLlowFreq=resample(GRLz,1,10);

%Creating NaN matrix with the lenght of the data 
newData=nan((((outmat(end,1)-outmat(1,1)))+1),15);
newData2=nan((((outmat(end,1)-outmat(1,1)))+1),15);

%Making frames from Python start at 1 
outmat(:,1)=outmat(:,1)-outmat(1,1)+1;

%find unique values for frame 
data=unique(outmat(:,1));

%finding unique colums
for i=1:length(data)
[datarows(i),~]=find(outmat(:,1)==data(i),1,'first');
end

%Creating a linear interpolate matrix from Python data 
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
 
NexusLlowFreq=NexusLlowFreq(1:length(newData));
NexusRlowFreq=NexusRlowFreq(1:length(newData));

%%
%Finding HS from Nexus at 100HZ and Interpolated Pyton data.
%Interpolate data is used to make sure that we dont take in consideration extras HS.
[LHSnexus,RHSnexus]= getEventsFromForces(NexusLlowFreq,NexusRlowFreq,100);
[LHSpyton,RHSpyton]= getEventsFromForces(newData(:,3),newData(:,2),100);

%% 
%localication of HS
locRHSpyton=find(RHSpyton==1);
locRHSnexus=find(RHSnexus==1);
locLHSpyton=find(LHSpyton==1);
locLHSnexus=find(LHSnexus==1);
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
GoodEvents=expData.data{j(w)}.adaptParams.Data(:,1); %Good events detected
GoodRHS=newData(locRindex,6);
GoodLHS=newData(locLindex,7); 
GoodRHS=GoodRHS(GoodEvents==1);
GoodLHS=GoodLHS(GoodEvents==1);
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
Rscale=newData(locLindex,10);
Lscale=newData(locLindex,11);
RscaleGood=Rscale(GoodEvents==1);
LscaleGood=Lscale(GoodEvents==1);
RtargetGood=(0.5./RscaleGood)*1000;
LtargetGood=(0.5./LscaleGood)*1000;    
end
%%
%plot of the alpha values. Tolerance indicade 
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
    %elseif GoodRHS(i)==1
        elseif abs(alphaRPytonGood(i)-RtargetGood(i))>=25
        plot(i,alphaRPytonGood(i),'k.','MarkerSize',15)
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
    %elseif GoodLHS(i)==1
     elseif abs(alphaLPytonGood(i)-LtargetGood(i))>=25
        plot(i,alphaLPytonGood(i),'k.','MarkerSize',15)
    end
  
end
title(['Alpha L leg' '(' ,subject ')' '  \alpha_(_n_)-\alpha^*<25' ])
axis tight
 end