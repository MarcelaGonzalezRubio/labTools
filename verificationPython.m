%SyncronizationPython 
% function a=verificationPython(subject,dynamic)
load(['Pyton1.mat'])
dynamic=0;
 
%Creating NaN matrix with the lenght of the data 
newData=nan((((outmat(end,1)-outmat(1,1)))+1),size(header,2));
newData2=nan((((outmat(end,1)-outmat(1,1)))+1),size(header,2));

%Making frames from Python start at 1 
outmat(:,1)=outmat(:,1)-outmat(1,1)+1;

%find unique values for frame 
data=unique(outmat(:,1));

%finding unique colums
for i=1:length(data)
[datarows(i),~]=find(outmat(:,1)==data(i),1,'first');
end
outmat=outmat(datarows,:);

%Calculating the gap's length for pyton data
gap=diff(outmat(:,1));
figure()
plot(gap,'b')

%Creating a linear interpolate matrix from Python data 
newData=interp1(outmat(:,1),outmat(:,1:end),[outmat(1,1):outmat(end,1)]);

%Creating a Matrix with NaN in gaps from Pyton 
for i=1:length(outmat)
  newData2(outmat(i,1),1:end)=outmat(i,:);
end


%%
%Finding HS Interpolated Pyton data.
%Interpolate data is used to make sure that we dont take in consideration extras HS.
[LHSpyton,RHSpyton]= getEventsFromForces(newData(:,3),newData(:,2),100);

%% 
%localication of HS
locRHSpyton=find(RHSpyton==1);
locLHSpyton=find(LHSpyton==1);
locRindex=find(newData2(:,4)==1);
locLindex=find(newData2(:,5)==1);
if length(locRindex)>length(locRHSpyton)
    display('Python detected more RHS')
% Delete extras HS deteted by Python 
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
else 
    display('Python did not detecte all HS')
end
%%
%Good strides
GoodRHS=newData2(locRindex,6);
GoodLHS=newData2(locLindex,7); 

%%
%find alpha value on time
alphaR_time=nan(length(newData2),1);
alphaL_time=nan(length(newData2),1);
alphaR_time(locRindex,1)=newData2(locRindex,8)*1000;
alphaL_time(locLindex,1)=newData2(locLindex,9)*1000;
%alpha values at HS 
alphaRPyton=newData2(locRindex,8)*1000;
alphaLPyton=newData2(locLindex,9)*1000;

% 
if dynamic ==1 
Rtarget=newData2(locRindex,10)*1000;
Ltarget=newData2(locLindex,11)*1000;
elseif dynamic== 0 %static target
Rscale=newData2(locRindex,10);
Lscale=newData2(locLindex,11);
Rtarget2=(0.25./Rscale)*1000;
Ltarget2=(0.25./Lscale)*1000;    
Rtarget=newData2(locRindex,16)*1000;
Ltarget=newData2(locLindex,17)*1000;
end
% LHip=newData2(locRindex,13)*1000;
% LAnk=newData2(locLindex,15)*1000;
hola=labTimeSeries(newData2,0,0.01,{'FrameNumber','Rfz','Lfz','RHS','LHS','RGORB','LGORB','Ralpha','Lalpha','Rscale','Lscale','RHIPY','LHIPY','RANKY','LANKY','targetR','targetL'});
%%
%plot of the alpha values. Tolerance indicade 
ystdRU=25*ones([length(GoodRHS),1])+Rtarget;
ystdRL=-25*ones([length(GoodRHS),1])+Rtarget;

ystdLU=25*ones([length(GoodLHS),1])+Ltarget;
ystdLL=-25*ones([length(GoodLHS),1])+Ltarget;
% ystdR2=25./Rscale;
% ystdL2=25./Lscale;
% end 
figure()
hold on 
% plot(1:length(GoodRHS),ystdRU,'--r',1:length(GoodRHS),ystdRL,'--r',1:length(GoodRHS),Rtarget,'r')
% errorbar(Rtarget2,ystdR,'rx');
% errorbar(Rtarget2,ystdR2,'rx');
for i=1:length(GoodRHS)
%         if abs((Rscale(i)*alphaRPyton(i))-500)<25
       if GoodRHS(i)==1
       %if abs(alphaRPyton(i)-Rtarget(i))<25
         plot(i,alphaRPyton(i),'o','MarkerSize',8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','g');
         GoodR(i,1)=1;
%         elseif abs((Rscale(i)*alphaRPyton(i))-500)>=25
       elseif GoodRHS(i)==0
%        elseif abs(alphaRPyton(i)-Rtarget(i))>=25
        plot(i,alphaRPyton(i),'o','MarkerSize',8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0]);
         GoodR(i,1)=0;
        end
    
end
 title(['Alpha R leg  \alpha_(_n_)-\alpha^*<25' ])
% title(['Alpha R leg \alpha_(_n_)-\alpha^*<25' ])
axis tight
%axis([0 1055 100 250])
figure()
hold on 

% errorbar(Ltarget2,ystdL,'rx');
% errorbar(Ltarget2,ystdL2,'rx');
% plot(1:length(GoodLHS),ystdLU,'--r',1:length(GoodLHS),ystdLL,'--r',1:length(GoodLHS),Ltarget,'r')
for i=1:length(GoodLHS)
    if GoodLHS(i)==1
% if abs((Lscale(i)*alphaLPyton(i))-500)<25
       %if abs(alphaLPyton(i)-Ltarget(i))<25
       plot(i,alphaLPyton(i),'o','MarkerSize',8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','g');
%        GoodL(i,1)=1;
    elseif GoodLHS(i)==0
%      elseif abs(alphaLPyton(i)-Ltarget(i))>=25
%  elseif abs((Lscale(i)*alphaLPyton(i))-500)>=25
        plot(i,alphaLPyton(i),'o','MarkerSize',8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0]);
         GoodL(i,1)=0;
    end
  
end
title(['Alpha L leg \alpha_(_n_)-\alpha^*<25' ])
axis tight
%title(['Alpha L leg' '(' ,subject   ')' '  \alpha_(_n_)-\alpha^*<25' ])
%title(['Alpha L leg \alpha_(_n_)-\alpha^*<25' ])
%axis([0 1055 100 250])
%