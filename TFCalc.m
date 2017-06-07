clear all
close all

Path = 'E:\FYP\TestData\' ;
file1 = '100avgs2cycles15ampchanel1TF3' ;
%file1 = '100avgschannel1TF' ;
file2 = 'IndividualChannels\1\1_channel';

%% Calculate transform function of system
load(strcat(Path, file1), 'avgIdealData', 'nSamplingFrequency', 'y', 'nAmplitude', 'timevec', 'nFrequency') 
load(strcat(Path,file2), 'outSignal', 'samplingRate', 'symbol');

[i, szout] = size(outSignal);

y=y.*nAmplitude;
szy = size(y);

signal = avgIdealData(:,1);
%start 5816 for 5 cycles, 5773 for 2 cycles
echostart = 5773;

echo = signal(echostart:echostart+szy(2)-1);

fft_x = (nSamplingFrequency/2700000)*[0:2700000-1];
time_y = (1/nSamplingFrequency)*[0:szy(2)-1];

time_echo = timevec(echostart:echostart+szy(2)-1) - timevec(echostart);

%time_echo = (1/nSamplingFrequency)*[0:szy(2)-1];

figure
plot(time_echo, echo)
hold on
plot(time_y, y)
xlabel('time [s]')
ylabel('Amplitude')
title('Waveform comparison')
legend('Echo', 'Ideal Pulse')
ylim([-1.5 1.5])

figure
plot(time_echo, echo./max(echo))
hold on
plot(time_y, y./max(y))
xlabel('time [s]')
ylabel('Amplitude')
title('Waveform comparison')
legend('Echo', 'Ideal Pulse')

echo = transpose(echo);

%-20 for 2 cycles, -50 for 5
echo = padarray(echo, [0 (2700000/2)-20]);
y = padarray(y, [0 (2700000/2)-20]);

echo_ft = (fft(echo));
y_ft = fft(y); 

echo_ft_n = abs(echo_ft)./max(abs(echo_ft));
y_ft_n = abs(y_ft)./max(abs(y_ft));

figure
plot(fft_x, mag2db(abs(echo_ft))-max(mag2db(abs(echo_ft))))
hold on
plot(fft_x, mag2db(abs(y_ft))-max(mag2db(abs(y_ft))))
xlim([0 10e6]);
xlabel('Frequency [Hz]')
ylabel('Amplitude [dB]') 
legend('Probe/Medium response','Ideal Signal')
title('Power Spectrums')

figure
plot(fft_x, mag2db(abs(echo_ft))./max(mag2db(abs(echo_ft))))
hold on
plot(fft_x, mag2db(abs(y_ft))./max(mag2db(abs(y_ft))))
xlim([0 10e6]);
xlabel('Frequency [Hz]')
ylabel('Amplitude [dB]') 
legend('Probe/Medium response','Ideal Signal')
title('Power Spectrums')

figure
plot(fft_x, mag2db(abs(echo_ft)))
hold on
plot(fft_x, mag2db(abs(y_ft)))
xlim([0 10e6]);
xlabel('Frequency [Hz]')
ylabel('Amplitude [dB]') 
legend('Probe/Medium response','Ideal Signal')
title('Power Spectrums')

figure
plot(fft_x, (abs(echo_ft)))
hold on
% yyaxis right
plot(fft_x, (abs(y_ft)))
xlim([0 10e6]);
xlabel('Frequency [Hz]')
ylabel('Amplitude [dB]') 
legend('Probe/Medium response','Ideal Signal')
title('Power Spectrums')

TF = echo_ft./y_ft;
% TF = y_ft./echo_ft;
TF_x = (nSamplingFrequency/2700000)*[0:2700000-1];

figure
plot(TF_x, mag2db(abs(TF)))
title('System Transfer Function')
xlabel('Frequency [Hz]')
ylabel('Amplitude [dB]')
xlim([0 10e6])

%% Apply TF to coded excitation

outSignal_ft = fft(outSignal);
outSignal_ft_x = (samplingRate/szout)*[0:szout-1];

figure
plot(outSignal_ft_x, abs(outSignal_ft))
hold on
yyaxis right
plot(TF_x, abs(TF))
xlim([0 samplingRate/2])
ylim([0 0.2])
title('Power spectrum of out signal')

%Downscale TF to match outsignal

TF=TF(1:270000);
TF=[TF, fliplr(TF)];


outSignal_TF_ft = TF.*outSignal_ft;

figure
plot(outSignal_ft_x, abs(outSignal_ft))
hold on
yyaxis right
plot(outSignal_ft_x, abs(outSignal_TF_ft)./max(abs(outSignal_TF_ft)))
xlim([0 samplingRate/2])
title('Power spectrums before and after applying TF')

outSignal_TF = real(ifft(outSignal_TF_ft)); 

idealecho = xcorr(outSignal,outSignal) ;
expectedecho = xcorr(outSignal,outSignal_TF) ;

idealecho = idealecho(539880:539905);
expectedecho = expectedecho(539880:539905);

time = 1/20 * [0:25];

figure
plot(time,idealecho)
ylim([-45000 45000])
xlim([0 1.2500])
ylabel('Amplitude')
xlabel('Time [\mus]')
title('Ideal Wave Form')

figure
plot(time,expectedecho)
ylim([-450 450])
xlim([0 1.2500])
ylabel('Amplitude')
xlabel('Time [\mus]')
%title('Expected Wave Form - Accounting for Probe Transfer Function')


for i = 1:540000
    if outSignal_TF(i) > 1.5e-3    
        outSignal_TF(i) = 1;
    elseif outSignal_TF(i) < -1.5e-3
        outSignal_TF(i) = -1;
    else
        outSignal_TF(i) = 0;
    end
end

expectedecho2 = xcorr(outSignal,outSignal_TF) ;
expectedecho2 = expectedecho2(539880:539905);

figure
plot(time,expectedecho2)
ylim([-45000 45000])
xlim([0 1.2500])
ylabel('Amplitude')
xlabel('Time [\mus]') 
%title('Expected Wave Form - Accounting for Binary Quantisation and Probe Transfer Function')

figure
plot(time,idealecho)
hold on
plot(time,expectedecho)
plot(time,expectedecho2)
legend('Sequence Autocorrelation', 'TF Applied', 'TF & Binary Quantisation Applied')
ylim([-7e4 7e4])
xlim([0 1.25])

figure
plot(time,idealecho/max(idealecho))
hold on
plot(time,expectedecho/max(expectedecho))
plot(time,expectedecho2/max(expectedecho2))
legend('Sequence Autocorrelation', 'TF Applied', 'TF & Binary Quantisation Applied')
title('Norm')
ylim([-1.75 1.75])
xlim([0 1.25])

%% Real Echo

path2 = 'G:\vary_voltage\repeat_' ;
load(strcat(path2, num2str(1),'\', num2str(4),'\', num2str(15000), '_ints' ),'xc','x_max') 
signal2 = xc( round(end/2) : round(end/2)+x_max);

%pad expected echo2 to show some noise

expectedecho2 = padarray(expectedecho2,[0 10]);

realecho = signal2(1+2272-size(expectedecho2,2)/2:2272+size(expectedecho2,2)/2);

figure
plot((0:x_max) /samplingRate *1e6,signal2)
ylim([-6e4 6e4])
xlim([0 140])
ylabel('Amplitude') 
xlabel('Time [\mus]')
%title('Location of Experimentally Obtained Echo')

time2 = 1/20 * [0:size(expectedecho2,2)-1];

figure
plot(time2,realecho)
xlim([0 2.2500])
%title('An Experimentally Obtained Echo')
ylabel('Amplitude')
xlabel('Time [\mus]')

figure
plot(time2,realecho./max(realecho))
hold on
plot(time2,expectedecho2./max(expectedecho2))
xlim([0 2.2500])
ylim([-1.5 1.5])
ylabel('Amplitude')
xlabel('Time [\mus]')
%title('Comparison of Experimentally Obtained and Theoretical Echoes')
legend('Experimentally Obtained Echo', 'Theoretically Predicted Echo')


