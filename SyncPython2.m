function results=SyncPython(subject)
%typeBiofeedback dymamics=1, Statics=0.
% subject='PDT06';
typeBiofeedback=1;
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
Stepsnexus=[];




results.locRindex=[];
results.locLindex=[];
results.locLindex2=[];
results.alphaRPyton=[];
results.alphaLPyton=[];
results.Rtarget=[];
results.Ltarget=[];
results.RscaleGood=[];
results.LscaleGood=[];
results.GoodRHS=[];
results.GoodLHS=[];
result.stepLengthAsym=[];

for p=1:length(condition)
    %  for p=2:2
    
    j=[];
    GRRz=[];
    GRLz=[];
    NexusRlowFreq=[];
    NexusLlowFreq=[];
    newData=[];
    newData2=[];
    datarows=[];
    locLHSnexus=[];
    locRHSnexus=[];
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
    alphaRPyton2=[];
    alphaLPyton2=[];
    alphaRPyton=[];
    alphaLPyton=[];
    Rtarget2=[];
    Ltarget2=[];
    Rtarget=[];
    Ltarget=[];
    Rscale=[];
    Lscale=[];
    RscaleGood=[];
    LscaleGood=[];
    GoodR=[];
    GoodL=[];
    Rtarget2Good=[];
    Ltarget2Good=[];
    alphaRnexus=[];
    alphaLnexusTemp=[];
    
    if strcmp(condition{p},'Gradual adaptation') || strcmp(condition{p},'Re-adaptation')
        w=w+1;
        
        load(['Pyton' num2str(w) '.mat'])
        z=expData.metaData.getConditionIdxsFromName(condition{p});
        j=adaptData.metaData.trialsInCondition{z};
        
        %         %Force plate data from Nexus
        GRRz=expData.data{j}.GRFData.Data(:,9);
        GRLz=expData.data{j}.GRFData.Data(:,3);
        
        %         %Converting force plate data from 960Hz to 120Hz
        NexusRlowFreq=resample(GRRz,1,9);
        NexusLlowFreq=resample(GRLz,1,9);
        
        %Creating NaN matrix with the lenght of the data
        newData=nan((((outmat(end,1)-outmat(1,1)))+1),2);
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
        %         gap=diff(outmat(:,1));
        %         figure()
        %         plot(gap,'b')
        
        %Creating a linear interpolate matrix from Pyton data
        newData=interp1(outmat(:,1),outmat(:,1:end),[outmat(1,1):outmat(end,1)]);
        
        %Creating a Matrix with NaN in gaps from Pyton
        for i=1:length(outmat)
            newData2(outmat(i,1),1:end)=outmat(i,:);
        end
        
        %%
%         LANK=expData.data{j}.markerData.Data(:,47);
        [RANK,time,auxLabel]=getDataAsVector(expData.data{j}.markerData,'RANKy');
        
        %Determination of crosscorrelation between Nexus at 100Hz and Interpolate
        %Pyton data
             figure()
            plot(RANK,'b')
            hold on
            plot(newData(:,16)*1000, 'r')
            
            [acor, lag]=xcorr(NexusRlowFreq,newData(:,2));
            [~,I]=max((acor));
            timeDiff=lag(I)
            
            if timeDiff<0
                newData=newData(abs(timeDiff)+1:end,1:end);
                newData2=newData2(abs(timeDiff)+1:end,1:end);
            elseif timeDiff>0
                newData=[zeros([timeDiff,size(header,2)]);newData];
                newData2=[zeros([timeDiff,size(header,2)]);newData2];
                
            end
            
            [acor, lag]=xcorr(RANK,newData(:,16));
            [~,I]=max((acor));
            timeDiff=lag(I)
            
            if timeDiff<0
                newData=newData(abs(timeDiff)+1:end,1:end);
                newData2=newData2(abs(timeDiff)+1:end,1:end);
            elseif timeDiff>0
                newData=[zeros([timeDiff,size(header,2)]);newData];
                newData2=[zeros([timeDiff,size(header,2)]);newData2];
                
            end
            figure()
            plot(RANK,'b')
            hold on
            plot(newData2(:,16)*1000, 'r')
        
        
        
        % Data=labTimeSeries(newData2,0,1/120,{'FrameNumber','Rfz','Lfz','RHS','LHS','RGORB','LGORB','Ralpha','Lalpha','R target','L target','RHIPY','LHIPY','RANKY','LANKY'});
        
        figure()
        plot(LANK,'b')
        hold on
        plot(newData2(:,15)*1000, 'r')
        %                 plot(newData2(:,16), 'g')
        %
        %         NexusLlowFreq=NexusLlowFreq(1:length(newData));
        %         NexusRlowFreq=NexusRlowFreq(1:length(newData));
        
        %%
        %Finding HS from Nexus at 100HZ and Interpolated Pyton data, interpolate
        %data is used to make sure that we dont take in consideration extras HS.
        [LHSnexus,RHSnexus,LTOnexus,RTOnexus]= getEventsFromForces(NexusLlowFreq,NexusRlowFreq,120);
        [LHSpyton,RHSpyton,LTOpyton,RTOpyton]= getEventsFromForces(newData(:,3),newData(:,2),120);
        
        %         if expData.data{j}.metaData.refLeg == 'L'
        %             locLHSnexus=adaptData.getParamInCond({'indSHSD'},condition{p});
        %             locRHSnexus=adaptData.getParamInCond({'indFHSD'},condition{p});
        %             locLTOnexus=adaptData.getParamInCond({'indSTOD'},condition{p});
        %             locRTOnexus=adaptData.getParamInCond({'indFTOD'},condition{p});
        %
        %         end
        %
        %         elseif expData.data{j}.metaData.refLeg == 'R'
        %             locLHSnexus=adaptData.getParamInCond({'indFHSD'},condition{p});
        %             locRHSnexus=adaptData.getParamInCond({'indSHSD'},condition{p});
        %         end
        %
        %                         figure()
        %                         plot(NexusRlowFreq,'b')
        %                         hold on
        %                         plot(newData(:,2), 'r')
        %                         plot(newData2(:,2), 'g')
        %         %                 plot(RTOnexus*100,'k')
        %         %                 plot(RTOpyton*100,'--m')
        %         %                 plot(newData2(:,6)*100,'--y')
        %         %                 legend('Nexus','Interpolate Pyton','Pyton','RHSnexus','RHSpyton','RHS python matrix')
        %                         title('Sync of R Force plate data')
        % %         %
        %                 figure()
        %                 plot(NexusLlowFreq,'b')
        %                 hold on
        % %                 plot(newData(:,3), 'r')
        %                 plot(newData2(:,3), 'g')
        % %                 plot(LTOnexus*100,'k')
        % %                 plot(LTOpyton*100,'--m')
        % %                 plot(newData2(:,7)*100,'--y')
        % %                 legend('Nexus','Interpolate Pyton','Pyton','LHSnexus','LHSpyton','LHS python matrix')
        %                 title('Sync of L Force plate data')
        %
        %% HEEL STRIKE
        %localication of HS==1);
        locLHSpyton=find(LHSpyton==1);
        locRHSpyton=find(RHSpyton==1);
        locRHSnexus=find(RHSnexus==1);
        locLHSnexus=find(LHSnexus==1);
        
        locRindex=find(newData2(:,4)==1);
        locLindex=find(newData2(:,5)==1);
        
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
            else
                locLindex(IsBad(1))=[];
            end
        end
        
        if length(locRHSnexus)<length(locRindex)
            FrameDiffR=[];
            IsBadR=[];
            while length(locRHSnexus)~=length(locRindex)
                diffLengthR=length(locRindex)-length(locRHSnexus);
                FrameDiffR=-locRindex(1:end-diffLengthR)+locRHSnexus;
                IsBadR=find(abs(FrameDiffR)>10);
                if isempty(IsBadR)
                    break
                else
                    locRindex(IsBadR(1))=[];
                end
            end
        end
        
        if length(locLHSnexus)<length(locLindex)
            FrameDiff=[];
            IsBad=[];
            while length(locLindex)~=length(locLHSnexus)
                diffLength=length(locLindex)-length(locLHSnexus);
                FrameDiff=-locLindex(1:end-diffLength)+locLHSnexus;
                IsBad=find(abs(FrameDiff)>10);
                if isempty(IsBad)
                    break
                else
                    locLindex(IsBad(1))=[];
                end
            end
        end
        
        if length(locRHSnexus)>length(locRindex)
            warning(['Gaps affected RHS detection  ' condition{p} ])
            
            while length(locRHSnexus)>length(locRindex)
                diffLengthR=-length(locRindex)+length(locRHSnexus);
                FrameDiffR=locRHSnexus(1:end-diffLengthR)-locRindex;
                IsBadR=find(FrameDiffR<=-10);
                if isempty(IsBadR)
                    break
                else
                    locfakeR=[locRindex(1:IsBadR-1);locRHSnexus(IsBadR(1));locRindex(IsBadR:end)];
                    locRindex=locfakeR;
                end
            end
        end
        if length(locLHSnexus)>length(locLindex)
            warning(['Gaps affected LHS detection  ' condition{p}])
            
            while length(locLHSnexus)>length(locLindex)
                diffLengthL=-length(locLindex)+length(locLHSnexus);
                FrameDiffL=locLHSnexus(1:end-diffLengthL)-locLindex;
                IsBadL=find(FrameDiffL<=-10);
                if isempty(IsBadL)
                    break
                else
                    locfakeL=[locLindex(1:IsBadL-1);locLHSnexus(IsBadL(1));locLindex(IsBadL:end)];
                    locLindex=locfakeL;
                end
            end
            
        end
        
        for i=1:length(locLindex)-1
            locLindex2(i,1)=locLindex(i+1);
        end
        %% TOE OFF
        
        %                 locLTOpyton=find(LTOpyton==1);
        %                 locRTOpyton=find(RTOpyton==1);
        %         %           locRTOnexus=find(RTOnexus==1);
        %         %                 locLTOnexus=find(LTOnexus==1);
        %                 locRindexTO=find(newData2(:,6)==1);
        %                 locLindexTO=find(newData2(:,7)==1);
        %
        %                 if length(locRindexTO)<length(locRTOpyton)
        %                     warning('No all the TO where detected')
        %                 end
        %
        %
        %                 %Delete extras HS deteted by Python
        %                 while length(locRTOpyton)~=length(locRindexTO)
        %                     diffLengthRTO=length(locRindexTO)-length(locRTOpyton);
        %                     FrameDiffRTO=locRindexTO(1:end-diffLengthRTO)-locRTOpyton;
        %                     IsBadRTO=find(FrameDiffRTO<=-10);
        %                     if isempty(IsBadRTO)
        %                         break
        %                     else
        %                         locRindexTO(IsBadRTO(1))=[];
        %                     end
        %                 end
        %
        %                 while length(locLTOpyton)~=length(locLindexTO)
        %                     diffLengthTO=length(locLindexTO)-length(locLTOpyton);
        %                     FrameDiffTO=locLindexTO(1:end-diffLengthTO)-locLTOpyton;
        %                     IsBadTO=find(FrameDiffTO<=-10);
        %                     if isempty(IsBadTO)
        %                         break
        %                     else
        %                         locLindexTO(IsBadTO(1))=[];
        %                     end
        %                 end
        %
        %                 if length(locRTOnexus)<length(locRindexTO)
        %                     FrameDiffRTO=[];
        %                     IsBadRTO=[];
        %                     while length(locRTOnexus)~=length(locRindexTO)
        %                         diffLengthRTO=length(locRindexTO)-length(locRTOnexus);
        %                         FrameDiffRTO=-locRindexTO(1:end-diffLengthRTO)+locRTOnexus;
        %                         IsBadRTO=find(abs(FrameDiffRTO)>10);
        %                         if isempty(IsBadRTO)
        %                             break
        %                         else
        %                             locRindexTO(IsBadRTO(1))=[];
        %                         end
        %                     end
        %                 end
        %
        %                 if length(locLTOnexus)<length(locLindexTO)
        %                     FrameDiffTO=[];
        %                     IsBadTO=[];
        %                     while length(locLindexTO)~=length(locLTOnexus)
        %                         diffLengthTO=length(locLindexTO)-length(locLTOnexus);
        %                         FrameDiffTO=-locLindexTO(1:end-diffLengthTO)+locLTOnexus;
        %                         IsBadTO=find(abs(FrameDiffTO)>10);
        %                         if isempty(IsBadTO)
        %                             break
        %                         else
        %                             locLindexTO(IsBadTO(1))=[];
        %                         end
        %                     end
        %                 end
        %
        %                 if length(locRTOnexus)>length(locRindexTO)
        %                     warning(['Gaps affected RHS detection  ' condition{p} ])
        %
        %                     while length(locRTOnexus)>length(locRindexTO)
        %                         diffLengthRTO=-length(locRindexTO)+length(locRTOnexus);
        %                         FrameDiffRTO=locRTOnexus(1:end-diffLengthRTO)-locRindexTO;
        %                         IsBadRTO=find(FrameDiffRTO<=-10);
        %                         if isempty(IsBadRTO)
        %                             break
        %                         else
        %                             locfakeRTO=[locRindexTO(1:IsBadRTO-1);locRTOnexus(IsBadRTO(1));locRindexTO(IsBadRTO:end)];
        %                             locRindexTO=locfakeRTO;
        %                         end
        %                     end
        %                 end
        %                 if length(locLTOnexus)>length(locLindexTO)
        %                     warning(['Gaps affected LHS detection  ' condition{p}])
        %
        %                     while length(locLTOnexus)>length(locLindexTO)
        %                         diffLengthLTO=-length(locLindexTO)+length(locLTOnexus);
        %                         FrameDiffLTO=locLTOnexus(1:end-diffLengthLTO)-locLindexTO;
        %                         IsBadLTO=find(FrameDiffLTO<=-10);
        %                         if isempty(IsBadLTO)
        %                             break
        %                         else
        %                             locfakeLTO=[locLindexTO(1:IsBadLTO-1);locLTOnexus(IsBadLTO(1));locLindexTO(IsBadLTO:end)];
        %                             locLindexTO=locfakeLTO;
        %                         end
        %                     end
        %
        %                 end
        %
        %                    for i=1:length(locRindexTO)-1
        %                     locRindexTO2(i,1)=locRindexTO(i+1);
        %                    end
        
        %%
        %Good strides
        GoodEvents=expData.data{j}.adaptParams.Data(:,1);
        locRindex=locRindex((GoodEvents)==1,1);
        locLindex=locLindex((GoodEvents)==1,1);
        locLindex2=locLindex2((GoodEvents)==1,1);
        %                 locRindexTO2=locRindexTO2((GoodEvents)==1,1);
        %                 locRindexTO=locRindexTO((GoodEvents)==1,1);
        %                 locLindexTO=locLindexTO((GoodEvents)==1,1);
        
        GoodRHS=newData2(locRindex,6);
        GoodLHS=newData2(locLindex,7);
        GoodLHS2=newData2(locLindex2,7);
        results.locLindex=[results.locLindex;locLindex];
        results.locRindex=[results.locRindex;locRindex];
        results.locLindex2=[results.locLindex2;locLindex2];
        results.GoodRHS=[results.GoodRHS;GoodRHS];
        results.GoodLHS=[results.GoodLHS;GoodLHS];
        
        
        %%
        %find alpha value on time
        stepLengthSlow=(newData2(locLindex2,9)-newData2(locLindex2,8))*1000;
        
        stepLengthFast=(newData2(locRindex,8)-newData2(locRindex,9))*1000;
        stepLengthDiff=stepLengthFast-stepLengthSlow;
        stepLengthAsym=stepLengthDiff./(stepLengthFast+stepLengthSlow);
        result.stepLengthAsym=[result.stepLengthAsym;stepLengthAsym];
        
        alphaR_time=nan(length(newData2),1);
        alphaL_time=nan(length(newData2),1);
        alphaR_time(locRindex,1)=newData2(locRindex,8)*1000;
        alphaL_time(locLindex,1)=newData2(locLindex,9)*1000;
        %alpha values at HS
        alphaRPyton=newData2(locRindex,8)*1000;
        alphaLPytonTemp=newData2(locLindex,9)*1000;
        alphaLPyton=newData2(locLindex2,9)*1000;
        %                 betaRpyton=newData2(locRindexTO2,10)*1000;
        %                 betaLpyton=newData2(locLindexTO,11)*1000;
        XRpyton=newData2(locRindex,9)*1000; % position of Left leg at RHS
        LRpyton=newData2(locLindex,8)*1000; % position of Rigth leg at LHS
        results.alphaRPyton=[results.alphaRPyton;alphaRPyton];
        results.alphaLPyton=[results.alphaLPyton;alphaLPyton];
        %
        if typeBiofeedback ==1
            Rtarget=newData2(locRindex,10)*1000;
            Ltarget=newData2(locLindex,11)*1000;
            Ltarget2=newData2(locLindex2,11)*1000;
            results.Rtarget=[results.Rtarget;Rtarget];
            results.Ltarget=[results.Ltarget;Ltarget];
            
            
        elseif typeBiofeedback== 0 %static target
%             Rscale=newData2(locRindex,10);
%             Lscale=newData2(locLindex,11);
%             
%             Rtarget2Good=(0.375./RscaleGood)*1000;
%             Ltarget2Good=(0.375./LscaleGood)*1000;
%             Rtarget=newData(locRindex,10)*1000;
%             Ltarget=newData(locLindex,11)*1000;
%             Ltarget2=newData(locLindex2,11)*1000;
%             results.Rtarget=[results.Rtarget;Rtarget];
%             results.Ltarget=[results.Ltarget;Ltarget];
%             results.RscaleGood=[results.RscaleGood;Rscale];
%             results.LscaleGood=[results.LscaleGood;Lscale];
        end
        
        
        %%
        %Comprobando si los pasos fueron clasificados de la manera correcta
        for i=1:length(GoodRHS)
            if abs(alphaRPyton(i)-Rtarget(i))<=37.5
                GoodR(i,1)=1;
            elseif abs(alphaRPyton(i)-Rtarget(i))>37.5
                GoodR(i,1)=0;
            elseif isnan(GoodRHS(i))
                GoodR(i,1)=NaN;
            end
            if GoodR(i)~=GoodRHS(i)
                display(['BAD LABEL RIGHT LEG ' num2str(i) ' STEP'])
            end
        end
        
        for i=1:length(GoodLHS)
            if abs(alphaLPytonTemp(i)-Ltarget(i))<=37.5
                GoodL(i,1)=1;
            elseif abs(alphaLPytonTemp(i)-Ltarget(i))>37.5
                GoodL(i,1)=0;
            elseif isnan(GoodLHS(i))
                GoodL(i,1)=NaN;
            end
            if GoodL(i)~=GoodLHS(i)
                display(['BAD LABEL LEFT LEG ' num2str(i) ' STEP'])
            end
        end
        
        
        %%
        alphaRnexus=adaptData.getParamInCond({'alphaFast'},condition{p});
        alphaLnexusTemp=adaptData.getParamInCond({'alphaTemp'},condition{p});
        alphaLnexus=adaptData.getParamInCond({'alphaSlow'},condition{p});
        
        %plot of the alpha values. Tolerance indicade
        ystdRU=37.5*ones([length(GoodRHS),1])+Rtarget;
        ystdRL=-37.5*ones([length(GoodRHS),1])+Rtarget;
        
        ystdLU=37.5*ones([length(GoodLHS),1])+Ltarget;
        ystdLL=-37.5*ones([length(GoodLHS),1])+Ltarget;
        
        ystdLU2=37.5*ones([length(GoodLHS2),1])+Ltarget2;
        ystdLL2=-37.5*ones([length(GoodLHS2),1])+Ltarget2;
        
        
        figure()
        hold on
        toleranceR=plot(1:length(GoodRHS),ystdRU,'--r',1:length(GoodRHS),ystdRL,'--r',1:length(GoodRHS),Rtarget,'r');
        
        
        GoodnexusR=NaN(length(GoodRHS),1);
        for i=1:length(GoodRHS)
            
            if GoodRHS(i)==1
                gPr=plot(i,alphaRPyton(i),'o','MarkerSize',8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','g');
                
                
            elseif GoodRHS(i)==0
                
                blackPr=plot(i,alphaRPyton(i),'o','MarkerSize',8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0]);
                
            end
            
            if abs(alphaRnexus(i)-Rtarget(i))<37.5
                RNr=plot(i,alphaRnexus(i),'o','MarkerSize',8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r');
                GoodnexusR(i,1)=1;
            elseif abs(alphaRnexus(i)-Rtarget(i))>=37.5
                blackNr=plot(i,alphaRnexus(i),'o','MarkerSize',8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b');
                GoodnexusR(i,1)=0;
            end
        end
        title('Alpha Fast leg' )
        gPr=gPr(1,1);
        toleranceR=toleranceR(1,1);
        blackPr=blackPr(1,1);
        RNr=RNr(1,1);
        blackNr=blackNr(1,1);
        legend([gPr,blackPr,RNr,blackNr,toleranceR],'Good Steps Python','Bad Steps Python','Good Steps Nexus','Bad Steps Nexus','Tolerance')
        axis tight
        
        figure()
        hold on
        toleranceL2=plot(1:length(GoodLHS2),ystdLU2,'--r',1:length(GoodLHS2),ystdLL2,'--r',1:length(GoodLHS2),Ltarget2,'r');
        GoodnexusL=NaN(length(GoodLHS2),1);
        for i=1:length(GoodLHS2)
            if GoodLHS2(i)==1
                
                gPL=plot(i,alphaLPyton(i),'o','MarkerSize',8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','g');
                
            elseif GoodLHS2(i)==0
                
                blackPL=plot(i,alphaLPyton(i),'o','MarkerSize',8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0]);
                
            end
            
            if abs(alphaLnexus(i)-Ltarget2(i))<37.5
                RNL= plot(i,alphaLnexus(i),'o','MarkerSize',8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r');
                GoodnexusL(i,1)=1;
            elseif abs(alphaLnexus(i)-Ltarget2(i))>=37.5
                
                blackNL= plot(i,alphaLnexus(i),'o','MarkerSize',8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b');
                GoodnexusL(i,1)=0;
            end
            %
        end
        title('Alpha Slow  leg' )
        gPL=gPL(1,1);
        blackPL=blackPL(1,1);
        RNL=RNL(1,1);
        blackNL=blackNL(1,1);
        toleranceL2=toleranceL2(1,1);
        legend([gPL,blackPL,RNL,blackNL,toleranceL2],'Good Steps Python','Bad Steps Python','Good Steps Nexus','Bad Steps Nexus','Tolerance')
        axis tight
        
        
        figure()
        hold on
        title('Alpha Slow Temp leg')
        toleranceL=plot(1:length(GoodLHS),ystdLU,'--r',1:length(GoodLHS),ystdLL,'--r',1:length(GoodLHS),Ltarget,'r');
        l=0;
        for i=1:length(GoodLHS)
            if GoodLHS(i)==1
                
                gPL=plot(i,alphaLPytonTemp(i),'o','MarkerSize',8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','g');
                
            elseif GoodLHS(i)==0
                
                blackPL=plot(i,alphaLPytonTemp(i),'o','MarkerSize',8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0]);
                GoodL(i,1)=0;
            end
            if abs(alphaLnexusTemp(i)-Ltarget(i))<37.5
                RNL= plot(i,alphaLnexusTemp(i),'o','MarkerSize',8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r');
            elseif abs(alphaLnexusTemp(i)-Ltarget(i))>=37.5
                blackNL= plot(i,alphaLnexusTemp(i),'o','MarkerSize',8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b');
                l=l+1;
            end
            
        end
        
        gPL=gPL(1,1);
        blackPL=blackPL(1,1);
        RNL=RNL(1,1);
        blackNL=blackNL(1,1);
        toleranceL=toleranceL(1,1);
        title(['Alpha Slow leg'])
        legend([gPL,blackPL,RNL,blackNL,toleranceL],'Good Steps Python','Bad Steps Python','Good Steps Nexus','Bad Steps Nexus','Tolerance')
        axis tight
        
    end
    %%
    
    numberSteps=adaptData.getParamInCond('Good',condition{p});
    %     stepLengthAsymnexus=adaptData.getParamInCond('stepLengthAsym',condition{p});
    
    StepsR2=NaN(length(numberSteps),1);
    StepsL2=NaN(length(numberSteps),1);
    Steps2=NaN(length(numberSteps),1);
    Steps3=NaN(length(numberSteps),1);
    
    if  strcmp(condition{p},'Gradual adaptation') || strcmp(condition{p},'Re-adaptation')
        StepsR2=GoodRHS;
        StepsL2=GoodLHS2;
        %         figure
        %         hold on
        for o=1:length(StepsR2)
            
            if  ((StepsR2(o)+StepsL2(o))/2)==1
                Steps2(o,1)=1;
                %                 blah=plot(o,stepLengthAsym(o),'o','MarkerSize',8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','g');
            else
                Steps2(o,1)=0;
                %                 blah2=plot(o,stepLengthAsym(o),'o','MarkerSize',8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','k');
            end
            %
            if  (( GoodnexusL(o)+ GoodnexusR(o))/2)==1
                Steps3(o,1)=1;
                %                         blah3=plot(o,stepLengthAsymnexus(o),'o','MarkerSize',8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r');
            else
                Steps3(o,1)=0;
                %                         blah2=plot(o,stepLengthAsymnexus(o),'o','MarkerSize',8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b');
            end
        end
        %         blah=blah(1,1);
        %         blah2=blah2(1,1);
        %         title('stepLengthAsym')
        %         legend([blah,blah2],'Good Steps','Bad Steps')
        %         axis tight
    end
    Steps=[Steps;Steps2];
    StepsR=[StepsR;StepsR2];
    StepsL=[StepsL;StepsL2];
    Stepsnexus=[Stepsnexus;Steps3];
    
end

%%
pData=adaptData.data;
labels={'TargetHitR', 'TargetHitL', 'TargetHit','TargetNexus'};
[aux,idx]=pData.isaLabel(labels);
if all(aux)
    adaptData.data.Data(:,idx)=[StepsR,StepsL,Steps,Stepsnexus];
else
    this=paramData([adaptData.data.Data,StepsR,StepsL,Steps,Stepsnexus],[adaptData.data.labels 'TargetHitR' 'TargetHitL' 'TargetHit' 'TargetNexus'],adaptData.data.indsInTrial,adaptData.data.trialTypes);
    adaptData=adaptationData(rawExpData.metaData,rawExpData.subData,this);
end
saveloc=[];
save([saveloc subject 'params.mat'],'adaptData');
save([saveloc subject 'Pyton2.mat'],'newData2','header','outmat','results','GoodLHS','GoodRHS','GoodnexusR','GoodnexusL')

end
