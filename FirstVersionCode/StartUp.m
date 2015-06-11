
% cd('/Users/carlysombric/Desktop/NewOGProcessing/newMethod/ paramFiles')
Smatrix=makeSMatrix;
subs=[{subFileList(Smatrix.CB)} {subFileList(Smatrix.SB)}];% {subFileList(Smatrix.BS)}];% {subFileList(Smatrix.YASV)}];
%subs=[{'OG11params.mat'} {'OG211params.mat'}]
params=[{'spatialContributionNorm2'} {'stepTimeContributionNorm2'} {'velocityContributionNorm2'} {'netContributionNorm2'}];


 close all

%conds=[{'TM base'} {'adaptation'} {'re-adaptation'}];
%conds=[{'TM base'}];
%conds=[{'TM base'} {'TM post'}];
%conds=[{'OG base'} {'OG post'}];
 %conds=[{'adaptation'}];
% conds=[{'Abrupt adaptation'} {'re-adaptation'}];
conds=[{'Abrupt adaptation'} {'re-adaptation'}];
%conds=[{'re-adaptation'}];


%adaptationData.plotAvgTimeCourse(subs,params,conds,5)% How to do the
%general timecourse plotting 
%adaptationData.plotAvgTimeCourse(subs,params,conds,20)% How to do the

%adaptationData.plotAvgTimeCourse(subs,params,conds,5, 1)% How to plot
%indicvidual subjects

 A=barGroups(Smatrix, params, {'CB','SB'});%, 'YA', 'YASV'}); %this is how the bar plots are made
 %close all

% Whole + Adaptation + Readaptation
%  conds=[];
 AvgTimeCourse_Whole(subs,params,conds,A) %This is the function that I use to do the rate calculations

% Cropped + Adaptation + Readaptation
%conds=[{'adaptation'} {'re-adaptation'}];
%AvgTimeCourse_Cropped(subs,params,conds,A) %This is the function that I use to do the rate calculations
