%SyncronizationPython 
% function a=verificationPython(subject,dynamic)
load(['Pyton1.mat'])
dynamic=0;
 
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
outmat=outmat(datarows,:);
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
locRindex=find(newData(:,4)==1);
locLindex=find(newData(:,5)==1);
if length(locRindex)<length(locRHSpyton)
locRindex=locRHSpyton;
locLindex=locLHSpyton;
end
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
GoodRHS=newData(locRindex+70,6);
GoodLHS=newData(locLindex+70,7); 

%%
%find alpha value on time
alphaR_time=nan(length(newData),1);
alphaL_time=nan(length(newData),1);
alphaR_time(locRindex,1)=newData(locRindex,8)*1000;
alphaL_time(locLindex,1)=newData(locLindex,9)*1000;
%alpha values at HS 
alphaRPyton=newData(locRindex,8)*1000;
alphaLPyton=newData(locLindex,9)*1000;

% 
if dynamic ==1 
Rtarget=newData(locRindex,10)*1000;
Ltarget=newData(locLindex,11)*1000;
elseif dynamic== 0 %static target
Rscale=newData(locRindex,10);
Lscale=newData(locLindex,11);
Rtarget=(0.5./Rscale)*1000;
Ltarget=(0.5./Lscale)*1000;    
end
%%
%plot of the alpha values. Tolerance indicade 
ystdR=25*ones([length(GoodRHS),1]);
ystdL=25*ones([length(GoodLHS),1]);
ystdR2=25./Rscale;
ystdL2=25./Lscale;
% end 
figure()
hold on 

errorbar(Rtarget,ystdR,'rx');
errorbar(Rtarget,ystdR2,'bx');

for i=1:length(GoodRHS)
%         if abs((Rscale(i)*alphaRPyton(i))-500)<25
       if GoodRHS(i)==1
        %if abs(alphaRPyton(i)-Rtarget(i))<25
         plot(i,alphaRPyton(i),'g.','MarkerSize',15)
         GoodR(i,1)=1;
%         elseif abs((Rscale(i)*alphaLPyton(i))-500)>=25
       elseif GoodRHS(i)==0
%        elseif abs(alphaRPyton(i)-Rtarget(i))>=25
        plot(i,alphaRPyton(i),'k.','MarkerSize',15)
         GoodR(i,1)=0;
        end
    
end
% title(['Alpha R leg' '(',subject ')' '   \alpha_(_n_)-\alpha^*<25' ])
title(['Alpha R leg \alpha_(_n_)-\alpha^*<25' ])
axis tight
figure()
hold on 
errorbar(Ltarget,ystdL,'rx');
errorbar(Ltarget,ystdL2,'bx');

for i=1:length(GoodLHS)
    if GoodLHS(i)==1
%  if abs((Lscale(i)*alphaLPyton(i))-500)<25
       %if abs(alphaLPyton(i)-Ltarget(i))<25
       plot(i,alphaLPyton(i),'g.','MarkerSize',15)
       GoodL(i,1)=1;
    elseif GoodLHS(i)==0
%      elseif abs(alphaLPyton(i)-Ltarget(i))>=25
%  elseif abs((Lscale(i)*alphaLPyton(i))-500)>=25
        plot(i,alphaLPyton(i),'k.','MarkerSize',15)
         GoodL(i,1)=0;
    end
  
end
% title(['Alpha L leg' '(' ,subject ')' '  \alpha_(_n_)-\alpha^*<25' ])
title(['Alpha L leg \alpha_(_n_)-\alpha^*<25' ])
axis tight
%