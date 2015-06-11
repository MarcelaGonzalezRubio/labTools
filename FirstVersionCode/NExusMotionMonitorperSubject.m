clear all 
load('BFP5.mat')
% 
% nexuscondition= {'pulse',...
%     'Rlimb',...
%     'Llimb',...
%     'RHIP',...
%     'LHIP',...
%     'RANK',...
%     'LANK',...
%     'Forcex',...
%     };

pulse=expData.data{2}.GRFData.getDataAsVector({['PulseFx']});

% angle=expData.data{2}.angleData.getDataAsVector({['RLimb'],['LLimb']});
% Rlimb=angle(:,1);
% Llimb=angle(:,2);
% 
% markers=expData.data{2}.markerData.getDataAsVector({['RHIPy'],['LHIPy'],['RANKy'],['LANKy']});
%  RHIP=markers(:,1);
%  LHIP=markers(:,2);
%  RANK=markers(:,3);
%  LANK=markers(:,4);
% 
% events=expData.data{2}.gaitEvents.getDataAsVector({['LHS'],['RHS'],['LTO'],['RTO']});
% SHS=events(:,1);
% FHS=events(:,2);
% STO=events(:,3);
% FTO=events(:,4);

Ntime = length(pulse);
Npulse = - pulse;
Npulse(Npulse < 2) = 0;


%  figure(1)
%  plot(Npulse);

[Npeaks,Nindex] = findpeaks(Npulse,'MinPeakDistance',158000);

% figure(2)
% plot(Ntime,Npulse);
% hold on
% plot(Nindex,Npeaks,'g');

load MM1

t2t = 0:(length(MM1)-1);
t2p = MM1(:,1);

% figure(3)
% plot(t2t,t2p);

[t2peak,t2index] = findpeaks(t2p,'MinPeakDistance',158000);

%compare t2index with Nindex(2), where these trials should start lining up
t2diff = Nindex(2)-t2index;
MM=[NaN((t2diff),9);MM1];

% figure(4)
% plot(Ntime,Npulse);
% hold on
% plot(t2t+t2diff,t2p,'g');

load MM2

t3t = 0:(length(MM2)-1);
t3p = MM2(:,1);

% figure(3)
% plot(t2t,t2p);

[t3peak,t3index] = findpeaks(t3p,'MinPeakDistance',158000);

%compare t2index with Nindex(2), where these trials should start lining up
t3diff = Nindex(4)-t3index;
MM=[MM;NaN((t3diff-length(MM)),9);MM2];
% figure(5)
% plot(Ntime,Npulse);
% hold on
% plot(t2t+t2diff,t2p,'g');
% plot(t3t+t3diff,t3p,'r');


load MM3

t4t = 0:(length(MM3)-1);
t4p = MM3(:,1);

% figure(3)
% plot(t2t,t2p);

[t4peak,t4index] = findpeaks(t4p,'MinPeakDistance',158000);

%compare t2index with Nindex(2), where these trials should start lining up
t4diff = Nindex(5)-t4index;
MM=[MM;NaN((t4diff-length(MM)),9);MM3];
% figure(5)
% plot(Ntime,Npulse);
% hold on
% plot(t2t+t2diff,t2p,'g');
% plot(t3t+t3diff,t3p,'r');
% plot(t4t+t4diff,t4p,'g');

load MM4

t5t = 0:(length(MM4)-1);
t5p = MM4(:,1);

% figure(3)
% plot(t2t,t2p);

[t5peak,t5index] = findpeaks(t5p,'MinPeakDistance',158000);

%compare t2index with Nindex(2), where these trials should start lining up
t5diff = Nindex(6)-t5index;
MM=[MM;NaN((t5diff-length(MM)),9);MM4];
% figure(5)
% plot(Ntime,Npulse);
% hold on
% plot(t2t+t2diff,t2p,'g');
% plot(t3t+t3diff,t3p,'r');
% plot(t4t+t4diff,t4p,'g');
% plot(t5t+t5diff,t5p,'r');


load MM5

t6t = 0:(length(MM5)-1);
t6p = MM5(:,1);

% figure(3)
% plot(t2t,t2p);

[t6peak,t6index] = findpeaks(t6p,'MinPeakDistance',158000);

%compare t2index with Nindex(2), where these trials should start lining up
t6diff = Nindex(7)-t6index;
MM=[MM;NaN((t6diff-length(MM)),9);MM5];
% figure(5)
% plot(Ntime,Npulse);
% hold on
% plot(t2t+t2diff,t2p,'g');
% plot(t3t+t3diff,t3p,'r');
% plot(t4t+t4diff,t4p,'g');
% plot(t5t+t5diff,t5p,'r');
% plot(t6t+t6diff,t6p,'g');

load MM6

t7t = 0:(length(MM6)-1);
t7p = MM6(:,1);

% figure(3)
% plot(t2t,t2p);

[t7peak,t7index] = findpeaks(t7p,'MinPeakDistance',158000);

%compare t2index with Nindex(2), where these trials should start lining up
t7diff = Nindex(8)-t7index;
MM=[MM;NaN(t7diff-length(MM),9);MM6];
MM=[MM;NaN((length(pulse)-length(MM)),9)];
figure(1)
% plot(Npulse);
% hold on
% plot(MM(:,1),'r');
% plot(t3t+t3diff,t3p,'r');
% plot(t4t+t4diff,t4p,'g');
% plot(t5t+t5diff,t5p,'r');
% plot(t6t+t6diff,t6p,'g');
% plot(t7t+t7diff,t7p,'r');
% 

