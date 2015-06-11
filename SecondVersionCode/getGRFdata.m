%This is how I am going to batch stuff

%% Contribution BATCHING

cd('F:\HARD DRIVE DULCE\Exp0002\ReprocessData\Params file')
Smatrix=makeSMatrix;
subs=[subFileList(Smatrix.CG)];% subFileList(Smatrix.BS) subFileList(Smatrix.D)];


h = waitbar(0,'Please wait...');
for i=1:numel(subs)

    cd(['F:\HARD DRIVE DULCE\Exp0002\ReprocessData\' num2str(subs{i}(1:end-10)) '\Session 1']) 

    load([num2str(subs{i}(1:end-10)) '.mat'])
    MetaData=expData.metaData;
    subject=expData.subData.ID;
    z=expData.metaData.getConditionIdxsFromName('re-adaptation');
    j=expData.metaData.trialsInCondition{z};
    j=j(1);
    Angle=expData.data{j}.angleData;
    GRF=expData.data{j}.GRFData;
    
    figure()
    subplot(2,1,1),plot(Angle.Data(:,2),'r')
    hold on
    subplot(2,1,1),plot(Angle.Data(:,1))
    labels=Angle.labels;
    legend(labels{2},labels{1})
    subplot(2,1,2),plot(GRF.Data(:,15))
    title(subject)
    
    save([subject 'GRF.mat'],'Angle','GRF','MetaData')
    display('one loop done')
    i
    waitbar(i/numel(subs))
end
cd('F:\HARD DRIVE DULCE\Exp0002\ReprocessData\Params file')
close(h)
