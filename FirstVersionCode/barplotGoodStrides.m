function barplotGoodStrides
load('PDT06Pyton2.mat')
load('PDT06info.mat')
% load('PST09Pyton2.mat')
% load('PST09info.mat')
Subject=info.ID
% BLeftGoodSteps=0;
% for i=1:1:length(GoodL)   
% if GoodL(i,1)==1
% BLeftGoodSteps=BLeftGoodSteps+1;
% else
% BLeftGoodSteps=BLeftGoodSteps;  
% end 
% end
% 
% BRigthGoodSteps=0;
% for i=1:1:length(GoodR)   
% if GoodR(i,1)==1
% BRigthGoodSteps=BRigthGoodSteps+1;
% else
% BRigthGoodSteps=BRigthGoodSteps;
% end 
% end
% 
% 
% BLeft=GoodL(:,1);
% BLeft(find(isnan(BLeft)))=[];
% BlegthLeft=length(BLeft);
% 
% BRigth=GoodR(:,1);
% BRigth(find(isnan(BRigth)))=[];
% BlegthRigth=length(BRigth);
% 
% L(1,1)=BLeftGoodSteps/BlegthLeft;
% R(1,1)=BRigthGoodSteps/BlegthRigth;

% ALeftGoodSteps=0;
% ARigthGoodSteps=0;
% for i=1:1:length(AdaptationL)   
% if AdaptationL(i,2)==1
% ALeftGoodSteps=ALeftGoodSteps+1;
% else
% ALeftGoodSteps=ALeftGoodSteps;
% end
% end
% for n=1:1:length(AdaptationR)
% if AdaptationR(n,2)==1
%  ARigthGoodSteps=ARigthGoodSteps+1;   
% else
%  ARigthGoodSteps=ARigthGoodSteps;   
% end
% end
% 
% ALeft=AdaptationL(:,2);
% ALeft(find(isnan(ALeft)))=[];
% AlegthLeft=length(ALeft);
% 
% ARigth=AdaptationR(:,2);
% ARigth(find(isnan(ARigth)))=[];
% AlegthRigth=length(ARigth);
% 
% L(2,1)=ALeftGoodSteps/AlegthLeft;
% R(2,1)=ARigthGoodSteps/AlegthRigth;
% 
% SLeftGoodSteps=0;
% SRigthGoodSteps=0;
% for i=1:1:length(SplitL)   
% if SplitL(i,2)==1
% SLeftGoodSteps=SLeftGoodSteps+1;
% else
% SLeftGoodSteps=SLeftGoodSteps;
% end
% end
% 
% for n=1:1:length(SplitR)
% if SplitR(n,2)==1
%     SRigthGoodSteps=SRigthGoodSteps+1;   
% else
%     SRigthGoodSteps=SRigthGoodSteps;     
% end
% end
% 
% SLeft=SplitL(:,2);
% SLeft(find(isnan(SLeft)))=[];
% SlegthLeft=length(SLeft);
% 
% SRigth=SplitR(:,2);
% SRigth(find(isnan(SRigth)))=[];
% SlegthRigth=length(SRigth);
% 
% L(3,1)=SLeftGoodSteps/SlegthLeft;
% R(3,1)=SRigthGoodSteps/SlegthRigth;
% 
% RLeftGoodSteps=0;
% RRigthGoodSteps=0;
% for i=1:1:length(PostCatch_GoodL)   
% if PostCatch_GoodL(i,1)==1
%  RLeftGoodSteps=RLeftGoodSteps+1;
% else
%    RLeftGoodSteps=RLeftGoodSteps;
% end 
% end
% for i=1:1:length(PostCatch_GoodR)   
% if PostCatch_GoodR(i,1)==1
%    RRigthGoodSteps=RRigthGoodSteps+1;
% else
%    RRigthGoodSteps=RRigthGoodSteps;  
% end
% end
% 
% RLeft=PostCatch_GoodL(:,1);
% RLeft(find(isnan(RLeft)))=[];
% RlegthLeft=length(RLeft);
% 
% RRigth=PostCatch_GoodR(:,1);
% RRigth(find(isnan(RRigth)))=[];
% RlegthRigth=length(RRigth);
% 
% L(2,1)=RLeftGoodSteps/RlegthLeft;
% R(2,1)=RRigthGoodSteps/RlegthRigth;

Rnexus=sum(GoodnexusR==1)/length(GoodnexusR);
Lnexus=sum(GoodnexusL==1)/length(GoodnexusL);
Rpython=sum(GoodRHS==1)/length(GoodRHS);
Lpython=sum(GoodLHS==1)/length(GoodLHS);
a=(Rnexus-Rpython)*100
b=(Lnexus-Lpython)*100

errorR=(sum(GoodRHS==1)-sum(GoodnexusR==1))/sum(GoodnexusR==1)
errorL=(sum(GoodLHS==1)-sum(GoodnexusL==1))/sum(GoodnexusL==1)

%figure()
% for i=1:1
% hold on
% bar((1:1)+(.5+.5.*i),Rnexus,0.2,'FaceColor',[.8,.8,.8])
% bar((1:1)+(.7+.5*i),Lnexus,0.2,'FaceColor',[.0,.36,.6])
% end
% condition={'Gradual Adaptation' ,'Gradual Adaptation'};
% axis([1.8 2.5 0 1])
% xTickPos=2.1:.5:2*length(condition);
% set(gca,'XTick',xTickPos,'XTickLabel',condition)
% 
% legend('Fast Leg' ,'Slow Leg')
% title(['Good Steps' '(',Subject ') nexus'])
% figure()
% 
% for i=1:1:1
% hold on
% bar((1:1)+(.5+.5.*i),Rpython(i,1),0.2,'FaceColor',[.8,.8,.8])
% bar((1:1)+(.7+.5*i),Lpython(i,1),0.2,'FaceColor',[.0,.36,.6])
% end
% condition={'Gradual Adaptation','Re-adaptation'};
% legend( 'Fast Leg','Slow Leg')
% title(['Good Steps' '(',Subject ') python'])
% axis([1.8 2.5 0 1])
% xTickPos=2.1:.5:2*length(condition);
% set(gca,'XTick',xTickPos,'XTickLabel',condition)

end