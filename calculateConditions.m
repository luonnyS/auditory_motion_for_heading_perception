function calculateConditions()
% {visualDegree visualDistance visualTime, ...
% auditoryDegree auditoryDistance auditoryTime,...
% sourceNum sourceDegree(:) sourceDistance(:) sourceHead(:)}

global TRIALINFO
global AUDITORY
global VISUAL

% TRIALINFO.stimulusType = [0 1 2]; % 0 for visual only, 1 for auditory only, 2 for both provided

% 0 for segregation
inteCondition0 = {};
% index0 = {};
if any(ismember(TRIALINFO.intergration,0))
    audiCon = cell(length(AUDITORY.sourceNum),4);
    for i=1:length(AUDITORY.sourceNum)
        audiCon(i,:) = {AUDITORY.sourceNum{i},AUDITORY.sourceDegree(i),AUDITORY.sourceDistance(i),AUDITORY.sourceHeading{i}};
    end
    con0AtS = [sortrows(repmat(AUDITORY.headingTime',size(audiCon,1),1)),...
        repmat(audiCon,length(AUDITORY.headingTime),1)];
    con0AdAtS = [sortrows(repmat(AUDITORY.headingDistance',size(con0AtS,1),1)),...
        repmat(con0AtS,length(AUDITORY.headingDistance),1)];
    con0AS = [sortrows(repmat(AUDITORY.headingDegree',size(con0AdAtS,1),1)),...
        repmat(con0AdAtS,length(AUDITORY.headingDegree),1)];
    con0AM = [repmat(con0AS,length(AUDITORY.MotionCoherence),1),...
        sortrows(repmat(AUDITORY.MotionCoherence',size(con0AS,1),1))];
    con0AC = [repmat(con0AM,length(AUDITORY.coherence),1),...
        sortrows(repmat(AUDITORY.coherence',size(con0AM,1),1))];
    con0ACFI = [repmat(con0AC,length(AUDITORY.coherenceFrameInitial),1),...
        sortrows(repmat(AUDITORY.coherenceFrameInitial',size(con0AC,1),1))];
    con0ACFD = [repmat(con0ACFI,length(AUDITORY.coherenceFrameDuration),1),...
        sortrows(repmat(AUDITORY.coherenceFrameDuration',size(con0ACFI,1),1))];

    con0VdVt = [sortrows(repmat(VISUAL.headingDistance',length(VISUAL.headingTime),1)),...
        repmat(VISUAL.headingTime',length(VISUAL.headingDistance),1)];
    con0V = [sortrows(repmat(VISUAL.headingDegree',size(con0VdVt,1),1)),...
        repmat(con0VdVt,length(VISUAL.headingDegree),1)];
    
    % visual only
    con0T0 = cat(2,con0V,num2cell(nan(size(con0V,1),size(con0ACFD,2))));
    
    % auditory only
    con0T1 = cat(2,num2cell(nan(size(con0ACFD,1),size(con0V,2))),con0ACFD);
    
    % both provided
    con0T2={};
    % reserved %
    
    if any(ismember(TRIALINFO.stimulusType,0))
        inteCondition0 = cat(1,inteCondition0,con0T0);
    end
    if any(ismember(TRIALINFO.stimulusType,1))
        inteCondition0 = cat(1,inteCondition0,con0T1);
    end
    if any(ismember(TRIALINFO.stimulusType,2))
        inteCondition0 = cat(1,inteCondition0,con0T2);
    end
end

% 1 for intergration
inteCondition1 = {};
% index1 = [];
if any(ismember(TRIALINFO.intergration,1))
    audiCon = cell(length(AUDITORY.sourceNum),4);
    for i=1:length(AUDITORY.sourceNum)
        audiCon(i,:) = {AUDITORY.sourceNum{i},AUDITORY.sourceDegree(i),AUDITORY.sourceDistance(i),AUDITORY.sourceHeading{i}};
    end
    con1AtS = [sortrows(repmat(TRIALINFO.headingTime',size(audiCon,1),1)),...
        repmat(audiCon,length(TRIALINFO.headingTime),1)];
    con1AdAtS = [sortrows(repmat(TRIALINFO.headingDistance',size(con1AtS,1),1)),...
        repmat(con1AtS,length(TRIALINFO.headingDistance),1)];
    con1AS = [sortrows(repmat(TRIALINFO.headingDegree',size(con1AdAtS,1),1)),...
        repmat(con1AdAtS,length(TRIALINFO.headingDegree),1)];
    con1AM = [repmat(con1AS,length(AUDITORY.MotionCoherence),1),...
        sortrows(repmat(AUDITORY.MotionCoherence',size(con1AS,1),1))];
    con1AC = [repmat(con1AM,length(AUDITORY.coherence),1),...
        sortrows(repmat(AUDITORY.coherence',size(con1AM,1),1))];
    con1ACFI = [repmat(con1AC,length(AUDITORY.coherenceFrameInitial),1),...
        sortrows(repmat(AUDITORY.coherenceFrameInitial',size(con1AC,1),1))];
    con1ACFD = [repmat(con1ACFI,length(AUDITORY.coherenceFrameDuration),1),...
        sortrows(repmat(AUDITORY.coherenceFrameDuration',size(con1ACFI,1),1))];

    
    con1VdVt = [sortrows(repmat(TRIALINFO.headingDistance',length(TRIALINFO.headingTime),1)),...
        repmat(TRIALINFO.headingTime',length(TRIALINFO.headingDistance),1)];
    con1V = [sortrows(repmat(TRIALINFO.headingDegree',size(con1VdVt,1),1)),...
        repmat(con1VdVt,length(TRIALINFO.headingDegree),1)];
    
    % visual only
    con1T0 = cat(2,con1V,num2cell(nan(size(con1V,1),size(con1ACFD,2))));
    
    % auditory only
    con1T1 = cat(2,num2cell(nan(size(con1ACFD,1),size(con1V,2))),con1ACFD);
    
    % both provided
    con1T2 = cat(2,con1ACFD(:,1:3),con1ACFD);
    
    if any(ismember(TRIALINFO.stimulusType,0))
        inteCondition1 = cat(1,inteCondition1,con1T0);
    end
    if any(ismember(TRIALINFO.stimulusType,1))
        inteCondition1 = cat(1,inteCondition1,con1T1);
    end
    if any(ismember(TRIALINFO.stimulusType,2))
        inteCondition1 = cat(1,inteCondition1,con1T2);
    end
end

TRIALINFO.trialConditions = [inteCondition0;inteCondition1];