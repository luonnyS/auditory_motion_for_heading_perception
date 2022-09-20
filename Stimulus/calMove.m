function [x,y,z,fx,fy,fz] = calMove(condition,refreshRate)
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