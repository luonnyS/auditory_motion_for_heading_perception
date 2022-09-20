function [y,z] = generateDimensionField(distance,degree,frustumLeft,frustumRight,frustumDepth)
if iscell(distance)
    distance = cell2mat(distance);
end
if iscell(degree)
    degree = cell2mat(degree);
end
posY = max(distance)*max(sind(degree))+ frustumRight;
negY = max(distance)*min(sind(degree)) + frustumLeft;
y = max(abs([posY,negY]))*2;
z = max(distance)*max(cosd(degree))+frustumDepth;