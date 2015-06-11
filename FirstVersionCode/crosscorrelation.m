
%Pulse for the differents trials
GRRz=expData.data{4}.GRFData.Data(:,9);
% GRRz2=expData.data{2}.GRFData.Data(:,9);
% GRRz3=expData.data{3}.GRFData.Data(:,9);

%Load MM data
FzR1=data(:,2);
% FzR2=MM2(:,2);
% FzR3=MM3(:,2);
%FzR4=MM4(:,2);

B=[];
% MbaselineSync=[]
 %for i=1:1:3

[acor(1), lag(1)]=xcorr(GRRz,FzR1);
[~,I(1)]=max(abs(acor(1)));
timeDiff=lag(I(1));
B=[B;NaN(timeDiff,1);FzR];
% MbaselineSync=[MbaselineSync;NaN((timeDiff(i))-lenght(MbaselineSync),10);MM(i)]; 

% figure()
% plot(GRRz,'b')
% hold on 
% plot(, '--r')

% %Trial 2
% [acor2, lag2]= xcorr(GRRz,FzR2);
% [~,d]=max(abs(acor2));
% timeDiff2=lag2(d);
% B=[B;NaN(timeDiff2-length(B),1);FzR2];
% MbaselineSync=[MbaselineSync;NaN((timeDiff2)-length(MbaselineSync),10);MM2];
% 

% 
% %Trial 3
% [acor3, lag3]= xcorr(GRRz,FzR3);
% [~,ds]=max(abs(acor3));
% timeDiff3=lag3(ds);
% B=[B;NaN(timeDiff3-length(B),1);FzR3];
% MbaselineSync=[MbaselineSync; NaN((timeDiff3)-length(MbaselineSync),10);MM3];
% end

figure()
plot(GRRz,'b')
hold on 
plot(B, '--r')
% figure()
% plot(GRRz3,'b')
% hold on 
% plot(j, '--r')

%Trial 4 
% [acor4, lag4]= xcorr(GRRz,FzR4);
% [~,ds4]=max(abs(acor4));
% timeDiff4=lag4(ds4);
% B=[B;NaN(timeDiff4,1);FzR4];

% figure()
% plot(GRRz,'b')
% hold on 
% plot(B, '--r')

%time vector for each trial
% t_base=expData.data{13}.gaitEvents.Time;

% %RHS and LHS 
% events_base=expData.data{13}.gaitEvents.getDataAsVector({['LHS'],['RHS']});
% RHS_base=[events_base(:,2)];
% LHS_base=[events_base(:,1)];

%creating matrix with time vectors 
% BaslineSync=labTimeSeries(MbaselineSync,0,0.001,{'Lz','Rz','targetR','targetL','Rgorb','Lgorb','Rhsp','Lhsp','RHS','LHS'});

%selecting the time for each RHS
% tBaseRHS=t_base(RHS_base==1);
% BaselineBFR=BaslineSync.getSample(tBaseRHS);

%selecting the time for each LHS
% tBaseLHS=t_base(LHS_base==1);
% BaselineBFL=BaslineSync.getSample(tBaseLHS);


% BaselineR=[BaselineBFR(:,2) BaselineBFR(:,3) BaselineBFR(:,5) BaselineBFR(:,7) BaselineBFR(:,9)];%BaselineBFR2(:,2) BaselineBFR2(:,3) BaselineBFR2(:,5) BaselineBFR2(:,7);BaselineBFR3(:,2) BaselineBFR3(:,3) BaselineBFR3(:,5) BaselineBFR3(:,7)];
% BaselineL=[BaselineBFL(:,1) BaselineBFL(:,4) BaselineBFL(:,6) BaselineBFL(:,8) BaselineBFL(:,10)];% BaselineBFL2(:,1) BaselineBFL2(:,4) BaselineBFL2(:,6) BaselineBFL2(:,8); BaselineBFL3(:,1) BaselineBFL3(:,4) BaselineBFL3(:,6) BaselineBFL3(:,8)];
% label={ 'Force' 'Target' 'Good' 'HSp' 'HS'} ;
% save('Syncronization.mat','BaselineR','BaselineL','label','Subject') 
