%--------------------------------------------------------------------------
%
% BIOMEDICAL SIGNAL PROCESSING
% JANUARY 20TH 2021 EXAM
% GIULIA PALLADINO
%
%--------------------------------------------------------------------------

clear
close all
clc

% File loading
load EEG.mat

eeg_adhd = EEG_ADHD;
eeg_sano = EEG_Control;

% Parmeters
fc = fs;                                            % sampling frequency
fny = fc/2;                                         % Nyquist frequency
t = 1/fc:1/fc:length(EEG_Control)/fc;               % time axis
nch = 7;                                            % number of channels

% Mean value detraction from the signals
for i=1:nch
    eeg_sano = eeg_sano - mean(eeg_sano);
    eeg_adhd = eeg_adhd - mean(eeg_adhd);
end

% EEG Theta Band: 3-7 Hz
Rp = 1;
Wp = 3/fny;

[b,a] = cheby1(6, Rp, Wp, 'high');
freqz(b,a,512,fc)

eeg_mtheta = filtfilt(b,a,eeg_adhd');
eeg_stheta = filtfilt(b,a,eeg_sano');

Rp = 1;
Wp = 7/fny;

[b,a] = cheby1(6, Rp, Wp, 'low');
freqz(b,a,512,fc)

eeg_mtheta = filtfilt(b,a,eeg_mtheta);
eeg_stheta = filtfilt(b,a,eeg_stheta);

eegmtetanorm = eeg_mtheta/std(eeg_adhd(:))/4;
eegstetanorm = eeg_stheta/std(eeg_sano(:))/4;

% EEG Beta Band: 14-30 Hz
Rp = 1;
Wp = 14/fny;

[b,a] = cheby1(6, Rp, Wp, 'high');
freqz(b,a,512,fc)

eeg_mbeta = filtfilt(b,a,eeg_adhd');
eeg_sbeta = filtfilt(b,a,eeg_sano');

Rp = 1;
Wp = 30/fny;

[b,a] = cheby1(6, Rp, Wp, 'low');
freqz(b,a,512,fc)

eeg_mbeta = filtfilt(b,a,eeg_mbeta);
eeg_sbeta = filtfilt(b,a,eeg_sbeta);

eegmbetanorm = eeg_mbeta/std(eeg_adhd(:))/4;
eegsbetanorm = eeg_sbeta/std(eeg_sano(:))/4;

% Signals' plot
% Control patient

figure(1)
subplot(3,1,1)
show_signal(eeg_sano, 'k', 'Original signal control patient', fc)
subplot(3,1,2)
show_signal(eegstetanorm', 'r', 'Theta rhythm control patient', fc)
subplot(3,1,3)
show_signal(eegsbetanorm', 'b', 'Beta rhythm control patient', fc)

% ADHD patient

figure(2)
subplot(3,1,1)
show_signal(eeg_adhd, 'k', 'Original signal ADHD patient', fc)
subplot(3,1,2)
Show_signal(eegmtetanorm', 'r', 'Theta rhythm ADHD patient', fc)
subplot(3,1,3)
show_signal(eegmbetanorm', 'b', 'Beta rhythm ADHD patient', fc)

% PSD and power ratio of the two bands (theta and beta)

for i=1:nch
   Pt = conv(eegmtetanorm(:,i).^2, ones(fc,1), 'same');
   Pteta(:,i) = Pt(fc/2:fc:end);
   Pb = conv(eegmbetanorm(:,i).^2, ones(fc,1), 'same');
   Pbeta(:,i) = Pb(fc/2:fc:end);
end

rapp_m = Pteta./Pbeta;

for i=1:nch
   Pt = conv(eegstetanorm(:,i).^2, ones(fc,1), 'same');
   Pteta(:,i) = Pt(fc/2:fc:end);
   Pb = conv(eegsbetanorm(:,i).^2, ones(fc,1), 'same');
   Pbeta(:,i) = Pb(fc/2:fc:end);
end

rapp_s = Pteta./Pbeta;

figure(3)
for i=1:nch
    ss = rapp_s(:,i)/max([rapp_m(:,i); rapp_s(:,i)]);
    sm = rapp_m(:,i)/max([rapp_m(:,i); rapp_s(:,i)]); 
    plot (ss+i,'g'), plot(sm+i, 'r')
    hold on
end

% Mean value and standard deviation of ratios for each signal
% Statistic test
for i=1:nch
   [H(i),p(i)] = ttest2(rapp_m(:,i),rapp_s(:,i));
   disp(['Theta-Beta Ratio - mean(std): ADHD ' num2str(mean(rapp_m(i,:))) '(' num2str(std(rapp_m(i,:))) ') - Control ' num2str(mean(rapp_s(i,:))) '(' num2str(std(rapp_s(i,:))) ') - t-test p ' num2str(p(i))])
end




