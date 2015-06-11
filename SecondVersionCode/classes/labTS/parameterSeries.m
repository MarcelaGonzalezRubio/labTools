classdef parameterSeries < labTimeSeries
    %parameterSeries  Extends labTimeSeries to hold adaptation parameters
    %   
    %parameterSeries properties:
    %   hiddenTime
    %   bad
    %   stridesTrial
    %   stridesInitTime
    %   description
    %
    %parameterSeries methods:
    %   idk
    %   
    
    properties
        hiddenTime        
    end
    properties(Dependent)
       bad
       stridesTrial
       stridesInitTime
       description
       trialTypes
    end
    properties(Hidden)
       description_={}; 
       trialTypes_={};
    end
    
    methods
        function this=parameterSeries(data,labels,times,description)            
            this@labTimeSeries(data,1,1,labels);
            this.hiddenTime=times;            
            if length(description)==length(labels)                
                this.description_=description; %Needs to be cell-array of same length as labels
            else
                error('paramtereSeries:constructor','Description input needs to be same length as labels')
            end       
        end
        
        function this=setTrialTypes(this,types)
            this.trialTypes_=types;
        end
       
        
        %% Getters for dependent variabls
        function vals=get.bad(this)
            if this.isaParameter('bad')
                vals=this.getDataAsVector('bad');
            elseif this.isaParameter('good')
                vals=this.getDataAsVector('good');
                vals=~vals;
            else
                %This should never be the case. Setting all values as good.
                vals=false(size(this.Data,1),1);
            end
        end
        function vals=get.stridesTrial(this)
            vals=this.getDataAsVector('trial');
        end
        function vals=get.stridesInitTime(this)
            vals=this.getDataAsVector('initTime');
        end
        function vals=get.description(this)
%            if isfield(this,'description_')
              vals=this.description_; 
%            else
%               vals=cell(size(this.labels)); 
%            end
        end
        function vals=get.trialTypes(this)  
%             if isfield(this,'trialTypes_')
               vals=this.trialTypes_;
%             else
%                 disp('trying to access trialTypes')
%                vals={}; 
%             end
        end
        
        
        %% I/O
        function [bool,idx]=isaParameter(this,labels) %Another name for isaLabel, backwards compatib
            [bool,idx]=this.isaLabel(labels);
        end
        
        function inds=indsInTrial(this,t)
            if nargin<2 || isempty(t)
                inds=[];
            else
                inds=cell(length(t),1);
                for ii=1:length(t)
                    inds{ii,1}=find(this.stridesTrial==t(ii));
                end
            end
        end
        
        function [data,auxLabel]=getParameter(this,label) %Backwards compat
            [data,~,auxLabel]=this.getDataAsVector(label);
        end                  
       
  
        %% Modifiers
        function newThis=cat(this,other)
            if size(this.Data,1)==size(other.Data,1)
                if isempty(this.description)
                    thisDescription=cell(size(this.labels));
                else
                    thisDescription=this.description;
                end
                if isempty(other.description)
                    otherDescription=cell(size(other.labels));
                else
                    otherDescription=other.description;
                end 
                newThis=parameterSeries([this.Data other.Data],[this.labels(:); other.labels(:)],this.hiddenTime,[thisDescription(:); otherDescription(:)]); 
            else
                error('parameterSeries:cat','Cannot concatenate series with different number of strides');
            end
        end
        
        function newThis=addStrides(this,other)
            %TODO: Check that the labels are actually the same
            if ~isempty(other.Data)
                aux=other.getDataAsVector(this.labels);
                if size(this.Data,2)==size(other.Data,2)                    
                    newThis=parameterSeries([this.Data; aux],this.labels(:),[this.hiddenTime; other.hiddenTime],this.description(:)); 
                else
                    error('parameterSeries:addStrides','Cannot concatenate series with different number of parameters.');
                end
            else
                newThis=this;
            end
        end
        
        function newThis=getDataAsPS(this,labels,strides)
            if nargin<2 || isempty(labels)
                labels=this.labels;
            end
            extendedLabels=[{'bad';'trial';'initTime'} ;labels(:)];
            extendedLabels=unique(extendedLabels); %To avoid repeating bad, trial, initTime
            [bool,idx]=this.isaLabel(extendedLabels);
            idx=idx(bool);
            if nargin<3 || isempty(strides)
               strides=1:size(this.Data,1); 
            end
            newThis=parameterSeries(this.Data(strides,idx),this.labels(idx),this.hiddenTime(strides),this.description(idx));
        end
        
        function newThis=appendData(this,newData,newLabels) %For back compat
            other=parameterSeries(newData,newLabels,this.hiddenTime,cell(size(newLabels)));
            newThis=cat(this,other);
        end
        
        %% Other functions that need redefining:
        function [F,f]=fourierTransform(~,~)
            error('parameterSeries:fourierTransform','You cannot do that!')
            F=[];
            f=[];
        end
        
        function newThis=resample(this) %the newTS is respected as much as possible, but forcing it to be a divisor of the total time range
            error('parameterSeries:resample','You cannot do that!')
            newThis=[];
        end
        
        function newThis=resampleN(this) %Same as resample function, but directly fixing the number of samples instead of TS
            error('parameterSeries:resampleN','You cannot do that!')
            newThis=[];
        end
        
        %% Display
        function [h,h1]=plotAlt(this,h,labels,plotHandles,color)
            if nargin<5
                color=[];
            end
            if nargin<4
                plotHandles=[];
            end
            if nargin<3
                labels=[];
            end
            if nargin<2
                h=[];
            end
            [h,h1]=this.plot(h,labels,plotHandles,[],color,1);
            ll=findobj(h,'Type','Line');
            set(ll,'LineStyle','None','Marker','.')
%             if nargin<2 || isempty(h)
%                 h=figure;
%             else
%                 figure(h)
%             end
%             N=length(this.labels);
%             if nargin<3 || isempty(labels)
%                 relData=this.Data;
%                 relLabels=this.labels;
%             else
%                [relData,~,relLabels]=this.getDataAsVector(labels); 
%                N=size(relData,2);
%             end
%             bad=this.bad;
%             for i=1:N
%                 h1(i)=subplot(ceil(N/2),2,i);
%                 T=1:length(bad);
%                 hold on
%                 plot(T(bad==0),relData(bad==0,i),'.')
%                 plot(T(bad==1),relData(bad==1,i),'x')
%                 ylabel(relLabels{i})
%                 hold off
%             end
             linkaxes(h1,'x')
                
        end
    end
    
end

