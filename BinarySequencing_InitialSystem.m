% HS5_PEsequence_v2.m
% written by Daniel Darvill 12/2016
% based on Tiepie LibTiepie examples and code by Frederic Cegla
% script sents and receives a Hanning windowed Toneburst with user
% controllable settings
%========================================================================
% 
for repeat = 1:11 

    clearvars -except repeat LibTiePie
    %display(repeat)
    
    bitlength = [2, 5, 8 ,9, 10 ,11, 12, 13, 14, 15, 16] ;
    
%% Acquistions settings
% (the only variables that have to be set by the user)
nRecLength = 64*96000;
nSamplingFrequency = 50e6;
nAmplitude = 1;
nFrequency = 5e6;
nCycles = 3;
nAverages = 1;
Sequence_bit=12;
TBextend=2;
pathname = 'C:\DanielD\TestData\ADC_Seq_A1_C1\';
filename ='_bit_ADCseqtest.mat';
%%

% Obtain Sequence waveform signals.
tbbitlength=nSamplingFrequency/nFrequency*nCycles;

%create A signals with Pseudo random sequence with receive intervals
%(Allows code to be open source and also hidden)
[signal1,OnOffSigA,OnOffSigB] = PE_sequence_Gen_Lic16(nFrequency, nCycles, Sequence_bit, tbbitlength, TBextend);

% prepare HS5 for sending data.
FC_get_HS5ready;  
 % loop through several acquistions   
    for count=1:1
       FC_HS5_getData;
       %pause(0.1)
    end
 % Disable output:
 gen.OutputOn = false; 
 % Close generator:
 clear gen;
% Close oscilloscope:
 clear scp;
 
 %%   
%     % Plot results:
    figure(1); % raw data
    plot(timevec, PE_seq_Data);
    %axis([0 (scp.RecordLength / scp.SampleFrequency) min(darRangeMin) max(darRangeMax)]);
    hold
    plot(timevec, Masked_data,'r');
    xlabel('Time [s]');
    ylabel('Amplitude [V]');
    title('Raw Data');
    legend('Raw','Masked');
    
    figure(2); %correlation output
    plot(timevec, Result,'r')
    title('Correlated');

    xlim([0 0.0001])

%save data
save(strcat(pathname, num2str(Sequence_bit), filename));

end
