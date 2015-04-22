
% cd('/Users/carlysombric/Desktop/NewOGProcessing/newMethod/ paramFiles')
Smatrix=makeSMatrix;
subs=[{subFileList(Smatrix.LN)}];% {subFileList(Smatrix.OASV)} {subFileList(Smatrix.YA)} {subFileList(Smatrix.YASV)}];
%subs=[{'OG11params.mat'} {'OG211params.mat'}]
params=[{'spatialContribution'} {'stepTimeContribution'} {'velocityContribution'} {'netContribution'}];


 close all

%conds=[{'TM base'} {'adaptation'} {'re-adaptation'}];
%conds=[{'TM base'}];
%conds=[{'TM base'} {'TM post'}];
%conds=[{'OG base'} {'OG post'}];
 %conds=[{'adaptation'}];
conds={'adaptation'};% {'re-adaptation'}]
%conds=[{'re-adaptation'}];


%adaptationData.plotAvgTimeCourse(subs,params,conds,5)% How to do the
%general timecourse plotting 
%adaptationData.plotAvgTimeCourse(subs,params,conds,20)% How to do the

%adaptationData.plotAvgTimeCourse(subs,params,conds,5, 1)% How to plot
%indicvidual subjects

 A=barGroupsFirstSteps(Smatrix, params, {'LN'}); %this is how the bar plots are made
 %close all

% Whole + Adaptation + Readaptation
% close conds=[{'adaptation'} {'adaptation'}];
 AvgTimeCourse_Whole(subs,params,conds,A) %This is the function that I use to do the rate calculations

% Cropped + Adaptation + Readaptation
%conds=[{'adaptation'} {'re-adaptation'}];
%AvgTimeCourse_Cropped(subs,params,conds,A) %This is the function that I use to do the rate calculations
