function show_signal(signal,str,tit,fsamp)

% Function created to plot multiple signals in one figure.
% Each signal is represented on a row in the final plot.

% Inputs:
% signal: signal to plot
% str: string for the signal's color in the plot
% tit: title of the plot
% fsamp: sampling frequency

[numChannels,len] = size(signal);           % number of channels to plot
t = 1/fsamp:1/fsamp:len/fsamp;              % define time axis

for i = 1:numChannels
    vv = signal(i,:);
    plot(t,vv+i,'Color',str)
    hold on
end

title(tit)
