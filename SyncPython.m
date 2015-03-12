% function SyncPython(subject,typeBiofeedback)
%typeBiofeedback dymamics=1, Statics=0.
subject='PST07';
typeBiofeedback=0;
load([subject 'params.mat'])
load([subject '.mat'])
load([subject 'RAW.mat'])

condition= adaptData.metaData.conditionName;
% This is new ~~~~~~~~~~~~~~~~~
%Needed in case subject did not perform one of the conditions
%in the condition list

condition=condition(find(~cellfun(@isempty,adaptData.metaData.trialsInCondition)));
w=0;
StepsR=[];
StepsL=[];
Steps=[];
% for p=1:length(condition)
    for p=1:2
    
    j=[];
    GRRz=[];
    GRLz=[];
    NexusRlowFreq=[];
    NexusLlowFreq=[];
    newData=[];
    newData2=[];
    datarows=[];
    LHSnexus=[];
    RHSnexus=[];
    LHSpyton=[];
    RHSpyton=[];
    locRHSpyton=[];
    locRHSnexus=[];
    locLHSpyton=[];
    locLHSnexus=[];
    locRindex=[];
    locLindex=[];
    IsBadR=[];
    IsBad=[];
    GoodEvents=[];
    GoodRHS=[];
    GoodLHS=[];
    data=[];
    alphaR_time=[];
    alphaL_time=[];
    alphaRPyton=[];
    alphaLPyton=[];
    alphaRPytonGood=[];
    alphaLPytonGood=[];
    Rtarget=[];
    Ltarget=[];
    RtargetGood=[];
    LtargetGood=[];
    Rscale=[];
    Lscale=[];
    RscaleGood=[];
    LscaleGood=[];
    GoodR=[];
    GoodL=[];
    Rtarget2Good=[];
    Ltarget2Good=[];
    
    if strcmp(condition{p},'Gradual adaptation') || strcmp(condition{p},'Re-adaptation')
%         w=w+1;
w=1;
        load(['Pyton' num2str(w) '.mat'])
        z=expData.metaData.getConditionIdxsFromName(condition{p});
        j=adaptData.metaData.trialsInCondition{z};
        
        %Force plate data from Nexus
        GRRz=expData.data{j}.GRFData.Data(:,9);
        GRLz=expData.data{j}.GRFData.Data(:,3);
        
        %Converting force plate data from 1000Hz to 100Hz
        NexusRlowFreq=resample(GRRz,1,10);
        NexusLlowFreq=resample(GRLz,1,10);
        
        %Creating NaN matrix with the lenght of the data
        newData=nan((((outmat(end,1)-outmat(1,1)))+1),size(header,2));
        newData2=nan((((outmat(end,1)-outmat(1,1)))+1),size(header,2));
        
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
        if timeDiff<0
            newData=newData(abs(timeDiff)+1:end,1:end);
            newData2=newData2(abs(timeDiff)+1:end,1:end);
        elseif timeDiff>0
            NexusRlowFreq =NexusRlowFreq(abs(timeDiff)+1:end,1:end);
            NexusLlowFreq =NexusLlowFreq(abs(timeDiff)+1:end,1:end);
        end
        
        hola=labTimeSeries(newData2,0,0.01,{'FrameNumber','Rfz','Lfz','RHS','LHS','RGORB','LGORB','Ralpha','Lalpha','Rscale','Lscale','RHIPY','LHIPY','RANKY','LANKY','targetR','targetL'});

%          figure()
%         plot(NexusRlowFreq,'b')
%         hold on
%         plot(newData(:,2), 'r')
%         plot(newData2(:,2), 'g')
%         legend('Nexus','Interpolate Pyton','Pyton')
%         title('Sync of R Force plate data')
%                 
%         figure()
%         plot(NexusLlowFreq,'b')
%         hold on
%         plot(newData(:,3), 'r')
%         plot(newData2(:,3), 'g')
%         legend('Nexus','Interpolate Pyton','Pyton')
%         title('Sync of L Force plate data')
        %
        NexusLlowFreq=NexusLlowFreq(1:length(newData));
        NexusRlowFreq=NexusRlowFreq(1:length(newData));
        
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
        
        locRindex=find(newData(:,4)==1);
        locLindex=find(newData(:,5)==1);
        
        if length(locRindex)<length(locRHSpyton)
            warning('No all the HS where detected')
        end
        
        
        %Delete extras HS deteted by Python
        while length(locRHSpyton)~=length(locRindex)
            diffLengthR=length(locRindex)-length(locRHSpyton);
            FrameDiffR=locRindex(1:end-diffLengthR)-locRHSpyton;
            IsBadR=find(FrameDiffR<=-10);
            if isempty(IsBadR)
                break
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
                break
                display('Something is WRONG!')
            else
                locLindex(IsBad(1))=[];
            end
        end
        
        %          display(['RHS nexus  ' num2str(length(locRHSpyton)) '   LHS nexus  ' num2str(length(locLHSpyton))])
        %          display(['RHS pyton  ' num2str(length(locRindex)) '   LHS pyton  ' num2str(length(locLindex))])
        if length(locRHSnexus)~=length(locRindex)
            warning(['Gaps affected RHS detection  ' condition{p} ])
            
             while length(locRHSnexus)~=length(locRindex)
                diffLengthR=-length(locRindex)+length(locRHSnexus);
                FrameDiffR=locRHSnexus(1:end-diffLengthR)-locRindex;
                IsBadR=find(FrameDiffR<=-10);
                if isempty(IsBadR)
                    break
                    display('Something is WRONG!')
                else
                    locfakeR=[locRindex(1:IsBadR-1);locRHSnexus(IsBadR(1));locRindex(IsBadR:end)];
                    locRindex=locfakeR;
                end
            end
        end
        if length(locLHSnexus)~=length(locLindex)
            warning(['Gaps affected LHS detection  ' condition{p}])
           
            while length(locLHSnexus)~=length(locLindex)
                diffLengthL=-length(locLindex)+length(locLHSnexus);
                FrameDiffL=locLHSnexus(1:end-diffLengthL)-locLindex;
                IsBadL=find(FrameDiffL<=-10);
                if isempty(IsBadL)
                    break
                    display('Something is WRONG!')
                else
                    locfakeL=[locLindex(1:IsBadL-1);locLHSnexus(IsBadL(1));locLindex(IsBadL:end)];
                    locLindex=locfakeL;
                end
            end
        
        end
        %
        %%
        %Good strides
        GoodEvents=expData.data{j}.adaptParams.Data(:,1);
        GoodRHS=newData2(locRindex,6);
        GoodLHS=newData2(locLindex,7);
        GoodRHS=GoodRHS(GoodEvents==1);
        GoodLHS=GoodLHS(GoodEvents==1);
        
        %%
        %find alpha value on time
        alphaR_time=nan(length(newData2),1);
        alphaL_time=nan(length(newData2),1);
        alphaR_time(locRindex,1)=newData2(locRindex,8)*1000;
        alphaL_time(locLindex,1)=newData2(locLindex,9)*1000;
        %alpha values at HS
        alphaRPyton=newData(locRindex,8)*1000;
        alphaLPyton=newData(locLindex,9)*1000;
        alphaRPytonGood=alphaRPyton(GoodEvents==1);
        alphaLPytonGood=alphaLPyton(GoodEvents==1);
        %
        if typeBiofeedback ==1
            Rtarget=newData2(locRindex,10)*1000;
            Ltarget=newData2(locLindex,11)*1000;
            RtargetGood=Rtarget(GoodEvents==1);
            LtargetGood=Ltarget(GoodEvents==1);
        elseif typeBiofeedback== 0 %static target
            Rscale=newData2(locRindex,10);
            Lscale=newData2(locLindex,11);
            RscaleGood=Rscale(GoodEvents==1);
            LscaleGood=Lscale(GoodEvents==1);
            Rtarget2Good=(0.25./RscaleGood)*1000;
            Ltarget2Good=(0.25./LscaleGood)*1000;
            Rtarget=newData(locRindex,16)*1000;
            Ltarget=newData(locLindex,17)*1000;
            RtargetGood=Rtarget(GoodEvents==1);
            LtargetGood=Ltarget(GoodEvents==1);
        end
       %% 
%Comprobando si los pasos fueron clasificados de la manera correcta
 if typeBiofeedback==1 
     for i=1:length(GoodRHS)
         if abs(alphaRPytonGood(i)-RtargetGood(i))<=25
             GoodR(i,1)=1;
         elseif abs(alphaRPytonGood(i)-RtargetGood(i))>25
             GoodR(i,1)=0;
         end
         if GoodR(i)~=GoodRHS(i)
             display(['BAD LABEL RIGHT LEG ' num2str(i) ' STEP'])
         end
     end
     
     for i=1:length(GoodLHS)
         if abs(alphaLPytonGood(i)-LtargetGood(i))<=25
             GoodL(i,1)=1;
         elseif abs(alphaLPytonGood(i)-LtargetGood(i))>25
             GoodL(i,1)=0;
         end
          if GoodL(i)~=GoodLHS(i)
             display(['BAD LABEL LEFT LEG ' num2str(i) ' STEP'])
         end
     end
 end
         
 if typeBiofeedback==0
     for i=1:length(GoodRHS)
         if abs(alphaRPytonGood(i)-RtargetGood(i))<=25
             GoodR(i,1)=1;
          elseif abs(alphaRPytonGood(i)-RtargetGood(i))>25
             GoodR(i,1)=0;
         end
         if GoodR(i)~=GoodRHS(i)
             display(['BAD LABEL RIGHT LEG ' num2str(i) ' STEP'])
         end
     end
     
      for i=1:length(GoodLHS)
          if abs(alphaLPytonGood(i)-LtargetGood(i))<=25
             GoodL(i,1)=1;
          elseif abs(alphaLPytonGood(i)-LtargetGood(i))>25
             GoodL(i,1)=0;
         end
         if GoodL(i)~=GoodLHS(i)
             display(['BAD LABEL LEFT LEG ' num2str(i) ' STEP'])
         end
     end
 end

    
    end
  %%  
    
    numberSteps=adaptData.getParamInCond('Good',condition{p});
    
    StepsR2=NaN(length(numberSteps),1);
    StepsL2=NaN(length(numberSteps),1);
    Steps2=NaN(length(numberSteps),1);
    
    if  strcmp(condition{p},'Gradual adaptation') || strcmp(condition{p},'Re-adaptation')
        StepsR2=GoodRHS;
        StepsL2=GoodLHS;
        for o=1:length(StepsR2)
            if  ((StepsR2(o)+StepsL2(o))/2)==1
                Steps2(o,1)=1;
            else
                Steps2(o,1)=0;
            end
        end
    end
    
    Steps=[Steps;Steps2];
    StepsR=[StepsR;StepsR2];
    StepsL=[StepsL;StepsL2];
    
end
%%
% pData=adaptData.data;
% labels={'TargetHitR', 'TargetHitL', 'TargetHit'};
% [aux,idx]=pData.isaLabel(labels);
% if all(aux)
%     adaptData.data.Data(:,idx)=[StepsR,StepsL,Steps];
% else
% this=paramData([adaptData.data.Data,StepsR,StepsL,Steps],[adaptData.data.labels 'TargetHitR' 'TargetHitL' 'TargetHit'],adaptData.data.indsInTrial,adaptData.data.trialTypes);
% adaptData=adaptationData(rawExpData.metaData,rawExpData.subData,this);
% end
% saveloc=[];
% save([saveloc subject 'params.mat'],'adaptData');
% end