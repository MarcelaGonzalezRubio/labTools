load('NBF15.mat')

%GRRz=expData.data{4}.GRFData.Data(:,9);
 GRRz=NexuslowFreq;
% GRRz2=expData.data{2}.GRFData.Data(:,9);
% GRRz3=expData.data{3}.GRFData.Data(:,9);
% GRR=[GRRz;GRRz2;GRRz3];
%load('MMall.mat')
FzR=newData(:,2);
% FzR2=MM2(:,2);
% FzR3=MM3(:,2);


% figure()
% plot(DGRRz,'b')
% hold on 
% plot(FzR, '--r')

%Trial 1
[acor, lag]=xcorr(GRRz,FzR);
[~,I]=max((acor));
timeDiff=lag(I);
 B=[];
 B=[FzR(abs(timeDiff):end)];

% [acor2, lag2]=xcorr(GRRz,B);
% [~,I2]=max(abs(acor2));
% timeDiff2=lag2(I2)

%Trial 2
% [acor2, lag2]= xcorr(GRR,FzR2);
% [~,d]=max(abs(acor2));
% timeDiff2=lag2(d);
% h=[NaN(timeDiff2,1);FzR2];
% B=[B;NaN(timeDiff2-length(B),1);FzR2];
%     
% %Trial 3
% [acor3, lag3]= xcorr(GRR,FzR3);
% [~,ds]=max(abs(acor3));
% timeDiff3=lag3(ds);
% j=[NaN(timeDiff3,1);FzR3];
% B=[B;NaN(timeDiff3-length(B),1);FzR3];

figure()
plot(GRRz,'b')
hold on 
plot(B, '--r')

% 
% t=GRRz(timeDiff2:end,1);
% R2=corrcoef([t,FzR2((1:length(t)),1)])
% 

% 
% t_base=expData.data{4}.gaitEvents.Time;
% events_base=expData.data{4}.gaitEvents.getDataAsVector({['LHS'],['RHS']});
% LHS_base=events_base(:,1);
% RHS_base=events_base(:,2);
% 
% MbaselineSync=[];
% MbaselineSync=[NaN((timeDiff),8);MM1]; 
% MbaselineSync=[MbaselineSync;NaN((timeDiff2-length(MbaselineSync)),8);MM2];
% MbaselineSync=[MbaselineSync;NaN((timeDiff3-length(MbaselineSync)),8);MM3];
% 
% tBaseRHS=t_base(RHS_base==1);
% BaslineSync=labTimeSeries(MbaselineSync,0,0.001,{'Lz','Rz','targetR','targetL','Rgorb','Lgorb','RHS','LHS'});
% BaselineBFR=BaslineSync.getSample(tBaseRHS);
% 
% tBaseLHS=t_base(LHS_base==1);
% BaselineBFL=BaslineSync.getSample(tBaseLHS);

% BaselineR=[BaselineBFR(:,2) BaselineBFR(:,4) BaselineBFR(:,6) BaselineBFR(:,7) BaselineBFR(:,10)];
% BaselineL=[BaselineBFL(:,3) BaselineBFL(:,5) BaselineBFL(:,8) BaselineBFL(:,9) BaselineBFL(:,11)];
% AdaptationR=[AdaptationBFR(:,2) AdaptationBFR(:,4) AdaptationBFR(:,6) AdaptationBFR(:,7) AdaptationBFR(:,10)];
% AdaptationL=[AdaptationBFL(:,3) AdaptationBFL(:,5) AdaptationBFL(:,8) AdaptationBFL(:,9) AdaptationBFL(:,11)];
% SplitR=[SplitBFR(:,2) SplitBFR(:,4) SplitBFR(:,6) SplitBFR(:,7) SplitBFR(:,10)];
% SplitL=[SplitBFL(:,3) SplitBFL(:,5) SplitBFL(:,8) SplitBFL(:,9) SplitBFL(:,11)];
% RadapR=[RadapBFR(:,2) RadapBFR(:,4) RadapBFR(:,6) RadapBFR(:,7) RadapBFR(:,10)];
% RadapL=[RadapBFL(:,3) RadapBFL(:,5) RadapBFL(:,8) RadapBFL(:,9) RadapBFL(:,11)];
% label={ 'Target' 'Good Step' 'Heel Strike Pos' 'Toe off Pos' 'Heel Strike'} ;
% save('Syncronization.mat','BaselineR','BaselineL','AdaptationR','AdaptationL','SplitR','SplitL','RadapR','RadapL','label','Subject') 