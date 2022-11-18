close all;
clear all;

FileName=('Z:\LQY Experiment\split\auditoryMotion_split3_2211181623.mat');
[pathstr,name]=fileparts(FileName);
load(fullfile(pathstr,name));
color=['b','k','r','p'];


%get number of parameters (X here) u need to fit in AUDITORY struct
%coh Analyzer: "coherence"; heading Analyzer: "headingDegree"
Parameter_need_to_fit='coherence';% "coherence" or "headingDegree"
%------------get modality-------
if ismember(0,TRIALINFO.stimulusType)
   Modality_need_to_fit=VISUAL;
   VisTrial=find(isnan(cell2mat(conditionIndex(:,4))));
   Vis_conditionIndex=conditionIndex(VisTrial,:);
   Vis_choice=choice(VisTrial,:);
end

if ismember(1,TRIALINFO.stimulusType)
   Modality_need_to_fit=AUDITORY;
   AudiTrial=find(isnan(cell2mat(conditionIndex(:,1))));
   Audi_conditionIndex=conditionIndex(AudiTrial,:);
   Audi_choice=choice(AudiTrial,:);
end
   
if ismember(2,TRIALINFO.stimulusType)
   Modality_need_to_fit=AUDITORY;
   CombTrial=find(isnan(cell2mat(conditionIndex(:,1)))==0 & isnan(cell2mat(conditionIndex(:,4)))==0);
   Comb_conditionIndex=conditionIndex(CombTrial,:);
   Comb_choice=choice(CombTrial,:);
end
% AudiTrial=find(isnan(cell2mat(conditionIndex(:,1))));
% Audi_conditionIndex=conditionIndex(AudiTrial,:);
% VisTrial=find(isnan(cell2mat(conditionIndex(:,4))));
% Vis_conditionIndex=conditionIndex(VisTrial,:);
% CombTrial=find(isnan(cell2mat(conditionIndex(:,1)))==0 & isnan(cell2mat(conditionIndex(:,4)))==0);
% Comb_conditionIndex=conditionIndex(CombTrial,:);
%-------------------------------
Flag_for_DDM=0;
modality=1;
switch Parameter_need_to_fit
    case 'coherence'
          Number_X=length(getfield(Modality_need_to_fit, Parameter_need_to_fit));
          X=sort(cell2mat(getfield(Modality_need_to_fit, Parameter_need_to_fit)));  
          for i=1:length(choice)
              if cell2mat(conditionIndex(i,4))<0
                 heading_degree(i,1)=1;
              else heading_degree(i,1)=2;
              end
          end
          choice_X_list{1}=[choice(:,1),cell2mat(conditionIndex(:,12))];
    case 'headingDegree'
          Number_X=length(getfield(Modality_need_to_fit, Parameter_need_to_fit));
          X=sort(cell2mat(getfield(Modality_need_to_fit, Parameter_need_to_fit))); 
          
          if ismember(0,TRIALINFO.stimulusType)
          choice_X_list{:,:,modality}=[Vis_choice(:,1),cell2mat(Vis_conditionIndex(:,1))];
          modality=modality+1;
          end
          if ismember(1,TRIALINFO.stimulusType)
          choice_X_list{:,:,modality}=[Audi_choice(:,1),cell2mat(Audi_conditionIndex(:,4))];
          modality=modality+1;
          end
         if ismember(2,TRIALINFO.stimulusType)
          choice_X_list{:,:,modality}=[Comb_choice(:,1),cell2mat(Comb_conditionIndex(:,1))];
         end
       
    otherwise disp('other value');
end

% a matrix to count rightchoice times in each X (heading discrimination) (choice_X_list{l}(i,1)==2)
% a matrix to count correct choice times in each X (coherence) (choice_X_list{l}(i,1)==2)

for l=1:size(choice_X_list,3)
    Counter=zeros(Number_X,1);
  for i=1:size(choice_X_list{l},1)
     for j=1:Number_X  
       if (choice_X_list{l}(i,2)==X(j)) &&  (choice_X_list{l}(i,1)==2)
          Counter(j)=Counter(j)+1;
       end
     end
  end


% fit psychometric curve       
aChoiceTimes(1:Number_X,:)=size(choice_X_list{l},1)/Number_X;
%aChoiceTimes=aChoiceTimes';
aPR = Counter./aChoiceTimes;
figureNum = 1;
aUniqueDeg=X';
fitData = [aUniqueDeg,aPR,aChoiceTimes];
    
[aBias,aThreshold] = cum_gaussfit_max1(fitData(1:end,:));
xi = min(aUniqueDeg):0.1:max(aUniqueDeg);
y_fit = cum_gaussfit([aBias,aThreshold],xi);
% plot   
if ishandle(figureNum); end; figure(figureNum); set(gcf,'color','white');
plot([0,0],[0,1],'-.k');
hold on
plot(aUniqueDeg,aPR,'*');
plot(xi,y_fit,'-','color',color(l));
set(gca, 'xlim',[min(aUniqueDeg),max(aUniqueDeg)],'ylim',[0 1])
xlabel('Coherence');
ylabel('Proportion of "right" choice');
title(['Participant ']);
text(6,0.3+l*0.15,sprintf('\\it\\mu_{psy} = \\rm%6.3g\\circ',aBias),'color',color(l))
text(6,0.25+l*0.15,sprintf('\\it\\sigma_{psy} = \\rm%6.3g\\circ', aThreshold),'color',color(l));
hold on
end
saveas(gcf,fullfile(pathstr,name),'jpg');
%csv file for DDM
if Flag_for_DDM==1
 for i=1:size(choice_X_list,1)
     if ((choice_X_list(i,2)<0.5) && (choice_X_list(i,1)==1)) || ((choice_X_list(i,2)>0.5) && (choice_X_list(i,1)==2))
        choice_X_list(i,3)=1;
     end

     if ((choice_X_list(i,2)<0.5) && (choice_X_list(i,1)==2)) || ((choice_X_list(i,2)>0.5) && (choice_X_list(i,1)==1))
        choice_X_list(i,3)=0;
     end
    
     if choice_X_list(i,2)==0.5
             seed=rand(1);
             if (seed<0.5 && choice_X_list(i,1)==1) || (seed>0.5 && choice_X_list(i,1)==2)
             choice_X_list(i,3)=1;
             else choice_X_list(i,3)=0;
             end
     end   
 end

 cohList=[choiceTime(:,1),choice_X_list(:,2),choice_X_list(:,3),choice_X_list(:,1)];
 colNames={'rt','coh','correct','trgchoice'};
 cohTable=array2table(cohList,'VariableNames',colNames);
 writetable(cohTable,'Z:\LQY Experiment\Cohtest\csvfile\auditoryMotion_CohtestC_2208221444.csv');
end