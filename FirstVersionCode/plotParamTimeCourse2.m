function plotParamTimeCourse2(this,label)   
            
            figureFullScreen
            
            figsz=[0 0 1 1];
            %in pixels:
%             vertpad = 30/scrsz(4); %padding on the top and bottom of figure
%             horpad = 60/scrsz(3);  %padding on the left and right of figure
%             
            % Set colors
            poster_colors;
            
            
            %find subplot size with width to hieght ratio of 4:1
            [rows,cols]=subplotSize(length(label),1,4);
            
            conds=find(~cellfun(@isempty,this.metaData.conditionName));
            nConds=length(conds);
            nPoints=size(this.data.Data,1);            
            rowind=1;
            colind=0;            
            for l=1:length(label)
                
                % Set colors order
                ColorOrder=[p_red; p_orange; p_fade_green; p_fade_blue; p_plum; p_green; p_blue; p_fade_red; p_lime; p_yellow].^(l^10);
                set(gca,'ColorOrder',ColorOrder);
                dataPoints=NaN(nPoints,nConds);
                for i=1:nConds
                    rawTrials=this.metaData.trialsInCondition{conds(i)};
                    if ~isempty(rawTrials)
                            %trials=find(ismember(cell2mat(this.metaData.trialsInCondition),rawTrials));
                            for t=rawTrials
                            inds=this.data.indsInTrial{t};
                            dataPoints(inds,i)=this.getParamInTrial(label(l),t);
                        end
                    end
                end
                %find graph location
                bottom=figsz(4)-(rowind*figsz(4)/rows)+vertpad;
                left=colind*(figsz(3))/cols+horpad;
                rowind=rowind+1;
                if rowind>rows
                    colind=colind+1;
                    rowind=1;
                end
                %subplot('Position',[left bottom (figsz(3)/cols)-2*horpad (figsz(4)/rows)-2*vertpad]);   
                hold on
                plot(dataPoints,'.','MarkerSize',15)  
                axis tight

                title([label{1},'-',label{2},' (',this.subData.ID ')']) 

            end
            condDes = this.metaData.conditionName;
            legend(condDes(conds)); %this is for the case when a condition number was skipped
        end
        
      