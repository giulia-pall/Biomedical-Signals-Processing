%--------------------------------------------------------------------------
%
% BIOMEDICAL SIGNALS PROCESSING
% JULY 15TH 2021 EXAM
% GIULIA PALLADINO
%
%--------------------------------------------------------------------------

clear
close all
clc

% File loading
head = readmatrix("head.txt");
pelvis = readmatrix("pelvis.txt");

acc_head = head(:,2:4);
acc_pelvis = pelvis(:,2:4);

fc = 100;

% Euclidean Norm
norm_head = vecnorm(acc_head');
norm_pelvis = vecnorm(acc_pelvis');

% Time axis
t = (0:length(acc_head)-1)'/fc;

% Signals plots
figure, subplot(211), plot(t, acc_head), hold on, plot(t, norm_head)
legend('x','y','z','Norm'), xlabel('time (s)'), ylabel('Acceleration (m/s^2)')
title('Head Acceleration')

subplot(212), plot(t, acc_pelvis), hold on, plot(t, norm_pelvis)
legend('x','y','z','Norm'), xlabel('time (s)'), ylabel('Acceleration (m/s^2)')
title('Pelvis Acceleration')

% Signal Filtering
wp = 20/(fc/2);         % Band pass frequency: up to 20 Hz, normalised
wt = 25/(fc/2);         % Attenuated band frequency: 25 Hz, normalised
rp = 1;                 % Band pass attenuation
rt = 15;                % Attenuated band attenuation

[n, wn] = buttord(wp, wt, rp, rt);
[b, a] = butter(n, wn, 'low');

figure, freqz(b, a, 10000, fc)

% Double filtering to avoid phase distorptions
norm_head_filt = filtfilt(b,a,norm_head);
norm_pelvis_filt = filtfilt(b,a,norm_pelvis);

% Filtered signal plot
figure, plot(t, norm_head_filt), hold on, plot(t, norm_pelvis_filt)
legend('norm head','norm pelvis'), xlabel('time (s)')
ylabel('Acceleration (m/s^2)'), title('Filtered Norms')

% Isolating a part of the signal containing 3 step's cycles
[~, loc] = findpeaks(norm_pelvis_filt, 'MinPeakHeight', 11, 'MinPeakDistance',0.4*fc);
int = loc(7):loc(13);
head2 = norm_head_filt(int);
pelvis2 = norm_pelvis_filt(int);

% Filtering out mean value 
head2 = head2-mean(head2);
pelvis2 = pelvis2-mean(pelvis2);

% PSD estimate of the input signal using Welch's overlapped segment
% averaging estimator.
[Ph, fh] = pwelch(head2, hamming(length(head2)), 0, length(head2)*10, fc);
[Pp, fp] = pwelch(pelvis2, hamming(length(pelvis2)), 0, length(pelvis2)*10, fc);

% PSD plot
figure, plot(fh, Ph), hold on, plot(fp, Pp)
legend('PSD head','PSD pelvis'), xlabel('Frequency (Hz)')
ylabel('PSD'), title('Signals'' PSD')

% Power ratio
rapp = sum(Ph)/sum(Pp);

disp('Power ratio:')
disp(rapp)

