clear all
close all
Data=zeros(100,100);
%for i=1:8
       y=zeros(100,100);
       %filei=['Noise' num2str(1000*(i-1)+100*(j-1)+l) '.wav'];   
       fs = 44.1e3;
       duration = 2;
       Octa=["1 octave"];
       %Roct=randperm(length(Octa));
       %oscillationFrequency = randi([100,900],1,1);% about amplitude of dB, recommend range 100 900
       oscillationFrequency = [20,50,100,200,300,400,500,600];
       %osiFreGroup=oscillationFrequency(randperm(length(oscillationFrequency)));
       %centerFreq = randi([1000,10000],1,1);% center frequency of
       %filter,recommend range 1000~20000
       centerFreq = [300,600,1200,2400,4800,9600];
       cenFreGroup=centerFreq(randperm(length(centerFreq)));
       cenFreGroup=cenFreGroup(1:length(Octa));
       %bw = Octa;% bandwidth of filter
      for i=1:length(centerFreq)
          for j=1:length(oscillationFrequency)
              octFilt = octaveFilter(centerFreq(i),Octa,'SampleRate',fs);
              y=pinknoise(duration*fs);
              y= octFilt(y);
              y=y.*reshape(sin((1:numel(y))*2*pi/(oscillationFrequency(j)*fs/1000)),size(y,1),size(y,2));
              plot(y);
              filei=['Z:\Gulab\noise generator\test\Noise' num2str((i-1)*length(oscillationFrequency)+j) '.wav'];
              audiowrite(filei,y,fs);
          end
      end
%       for i=1:8
%        octFilt = octaveFilter(cenFreGroup(i),bw(i),'SampleRate',fs);
% % %        Data={100*i+10*j+l, oscillationFrequency(j), centerFreq(l), Roct(i)}l
% % %        Data(l+,1)=100*i+10*j+l;
% % %        Data(i,2)=oscillationFrequency(j);
% % %        Data(i,3)=centerFreq(l);
% % %        Data(i,4)=Roct(i);
%        y = pinknoise(duration*fs);
%        y= octFilt(y);
%        y=y.*reshape(sin((1:numel(y))*2*pi/(osiFreGroup(i)*fs/1000)),size(y,1),size(y,2));
%        plot(y);
%        filei=['Z:\Gulab\noise generator\220806 20\Noise' num2str(i+16 ) '.wav'];
%        audiowrite(filei,y,fs);
%       end
%  ----------------------------average power spectral density-------------------
%   [~,freqVec,~,psd] = spectrogram(y,round(0.05*fs),[],[],fs);
%   meanPSD = mean(psd,2);
%   semilogx(freqVec,db(meanPSD,"power"))
%   xlabel('Frequency (Hz)')
%   ylabel('PSD (dB/Hz)')
%   title('Power Spectral Density of Pink Noise (Averaged)')
%   grid on

%end