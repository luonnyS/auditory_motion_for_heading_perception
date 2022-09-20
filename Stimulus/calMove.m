function [x,y,z,fx,fy,fz] = calMove(condition,refreshRate,Axis)
% Horizontal
if Axis==0
degree = condition(1);
distance = condition(2);
time = condition(3);

frameNum = time * refreshRate+1;

x = (1:frameNum)/frameNum*distance*sind(degree);
y = (1:frameNum)*0;
z = (1:frameNum)/frameNum*distance*-cosd(degree);

fx = zeros(1,frameNum);
fy = zeros(1,frameNum);
fz = -ones(1,frameNum);
end
% Vertical
if Axis==1
degree = condition(1);
distance = condition(2);
time = condition(3);

frameNum = time * refreshRate+1;

y = (1:frameNum)/frameNum*distance*sind(degree);
x = (1:frameNum)*0;
z = (1:frameNum)/frameNum*distance*-cosd(degree);

fy = zeros(1,frameNum);
fx = zeros(1,frameNum);
fz = -ones(1,frameNum);
end

