%This is how I am going to batch stuff

%% Contribution BATCHING

% close all
% clear all
% clc

%cd('C:\Users\dum5\Desktop\dulce\Exp0002\Subjects\Final matrix subjects')
%cd('C:\Users\dum5\Desktop\dulce\Exp0002\matrix')
% cd('E:\Exp0002\Subjects\Final matrix subjects')
cd('E:\Exp0002\Subjects\Biofeedback Pilots\Pilot static and dynamic targets params')

Smatrix=makeSMatrix;
subs=[subFileList(Smatrix.PST) subFileList(Smatrix.PDT)];% subFileList(Smatrix.D)];


h = waitbar(0,'Please wait...');
for i=1:numel(subs)
    %cd(['C:\Users\dum5\Desktop\dulce\Exp0002\Subjects\' num2str(subs{i}(1:end-10)) '\Session 1'])
    cd(['E:\Exp0002\Subjects\Biofeedback Pilots\' num2str(subs{i}(1:end-10))])
    %cd(['C:\Users\dum5\Desktop\dulce\Exp0002\YA\' num2str(subs{i}(1:end-10)) ])
    makeDataObject(subs{i}(1:end-10))
    %Feedback_test1_rev2(subs{i}(1:end-10))
    %Feedback_test1_rev2(subs{i}(1:end-10))
    display('one loop done')
    i
    waitbar(i/numel(subs))
end

close(h)
