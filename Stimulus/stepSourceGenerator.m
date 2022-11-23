function stepSourceGenerator(duration,steps,audioFreqMax,audioFreqMin,simulateType,initial,terminal,fileName,enableOscillation,oscillationFrequency,testEnable)
if nargin<1 || isempty(duration)
    duration =5; % second
end
if nargin<2 || isempty(steps)
    steps = 1;
end
if nargin <3 || isempty(audioFreqMax)
    audioFreqMax = 500;
end
if nargin <4 || isempty(audioFreqMin)
    audioFreqMin = 500; % will be disable in pure tone ( steps == 1)
end
if nargin<5 || isempty(simulateType)
    simulateType = 0; % 0 for two channel, 1 for left ear only, 2 for right ear only, 3 to simulate left to right by loudness, 4 to simulate right to left by loudness
end
if nargin<6 || isempty(initial)
    initial = 0.05;% second
end
if nargin<7 || isempty(terminal)
    terminal = 0.05; % second
end
if nargin<8 || isempty(fileName)
    fileName = 'stepTone500.wav';
elseif ~ischar(fileName)
    fileName = 'stepTone.wav';
    warning('Invalid file name, rename as ''stepTone.wav''');
end
if nargin<9 || isempty(enableOscillation)
    enableOscillation = true;
end
if nargin < 10 || isempty(oscillationFrequency) 
    if enableOscillation
        oscillationFrequency = 400; % default 100ms
    end
end
if nargin<10 || isempty(testEnable)
    testEnable = true;
end

if initial > terminal
    iniLonger = true;
else
    iniLonger = false;
end

sampleRate = 96000; % how many sample per second, 1000 to 384000 Hz, typically 8000, 11025, 22050, 44100, 48000, and 96000 Hz
% audioFrequency=linspace(audioFreqMin,audioFreqMax,sampleRate*duration);
frequencies =  linspace(audioFreqMin,audioFreqMax,steps);

% /: initial        _: mainbody      \: terminal
% ////____\\           ////____\\
%               ////____\\
% or
% //____\\\\           //____\\\\
%               //____\\\\
if iniLonger
    bodyDuration = (duration-terminal)/steps-initial;
else
    bodyDuration = (duration-initial)/steps-terminal;
end
if bodyDuration < 0
    error('terminal/initial too long');
end

y = zeros(round(sampleRate*duration),1); % matlab have some problam for big number, so I have to fun round() on it.
for i = 1:steps
    t = (1:round((bodyDuration+initial+terminal)*sampleRate))/sampleRate;
    tempY = sin(frequencies(i)*2*pi.*t)';
    initialAmp = linspace(0,1,sampleRate*initial);
    terminalAmp = linspace(1,0,sampleRate*terminal);
    tempY(1:length(initialAmp)) = tempY(1:length(initialAmp)).*initialAmp';
    tempY(end-length(terminalAmp)+1:end) = tempY(end-length(terminalAmp)+1:end).*terminalAmp';
    if i == 1
        y(1:round((bodyDuration+initial+terminal)*sampleRate)) = tempY;
    else
        if iniLonger
            y((i-1)*(bodyDuration+initial)*sampleRate+1 : (i*(bodyDuration+initial)+terminal)*sampleRate) = y((i-1)*(bodyDuration+initial)*sampleRate+1 : (i*(bodyDuration+initial)+terminal)*sampleRate) +tempY;
        else
            y((i-1)*(bodyDuration+terminal)*sampleRate+1 : (i*(bodyDuration+terminal)+initial)*sampleRate) = y((i-1)*(bodyDuration+terminal)*sampleRate+1 : (i*(bodyDuration+terminal)+initial)*sampleRate)+tempY;
        end
    end
end

if enableOscillation
    y=y.*reshape(sin((1:numel(y))*2*pi/(oscillationFrequency*sampleRate/1000)),size(y,1),size(y,2));
end
    
switch simulateType
    case 0
        yfin = y;
    case 1
        yfin = [y, zeros(size(y))];
    case 2
        yfin = [zeros(size(y)),y];
    case 3
        yleft = sin((1:round(duration*sampleRate))/round(duration*sampleRate)*pi/2);
        yright = cos((1:round(duration*sampleRate))/round(duration*sampleRate)*pi/2);
        yfin = [y.*yleft',y.*yright'];
    case 4
        yleft = cos((1:round(duration*sampleRate))/round(duration*sampleRate)*pi/2);
        yright = sin((1:round(duration*sampleRate))/round(duration*sampleRate)*pi/2);
        yfin = [y.*yleft',y.*yright'];
end
player = audioplayer(yfin, sampleRate);
if testEnable
    plot(yfin)
    player.play();
    pause(duration);
end
audiowrite(fileName, yfin, sampleRate);
end