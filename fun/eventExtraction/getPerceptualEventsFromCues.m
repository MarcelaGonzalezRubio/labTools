function [LHSstartCue, LHSstopCue, RHSstartCue, RHSstopCue] = getPerceptualEventsFromCues(datlog, infoLHSevent, infoRHSevent) 
 
% Grab auditory cues time from the datlog. This information has to be
% offset following the synchronization process between datlogs and Nexus
startCue = datlog.audioCues.start + datlog.dataLogTimeOffsetBest;
stopCue = datlog.audioCues.stop + datlog.dataLogTimeOffsetBest;

%% Compare the start and stop cue times to the events data to match the start and stop of perceptual trial 
LHSstartCue = zeros(length(infoLHSevent),1); 
RHSstartCue = zeros(length(infoRHSevent),1); 
LHSstopCue = zeros(length(infoLHSevent),1); 
RHSstopCue = zeros(length(infoRHSevent),1); 

    for t=1:length(startCue)
    
        [~,LidxI]=min(abs(infoLHSevent - startCue(t)));
        [~,LidxF]=min(abs(infoLHSevent - stopCue(t)));
        LHSstartCue(LidxI) = infoLHSevent(LidxI);
        LHSstopCue(LidxF) = infoLHSevent(LidxF);
    
        [~,RidxI]=min(abs(infoRHSevent - startCue(t)));
        [~,RidxF]=min(abs(infoRHSevent - stopCue(t)));
        RHSstartCue(RidxI) =  infoRHSevent(RidxI);
        RHSstopCue(RidxF) = infoRHSevent(RidxF);
    
    end

end

