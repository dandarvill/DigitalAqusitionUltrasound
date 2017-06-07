%**************************************************************************
%
% ab_io_std_sync.m               Created by Kostas Zarogoulidis, 10/01/2017
%                                Based on io_std_sync by Spectrum
%
%**************************************************************************
%
% This script performs A/B synchronisation between two M2i cards, with the 
% first being on standard replay mode, while the second on standard 
% recording mode. The Star-Hub option is used for the syncrhonisation.
%
% ANALOG OUTPUT (THE ONE SEND TO TRANS): outSignalsNoDump 
% 
% DIGITAL INPUT : DigDataExpanded 
%
% READ 'M2i.7020 Synchronous IO.docx' for detail explanation of
% parameters/variables
%**************************************************************************

% ints = [50 100 250 500 1000 2500 5000 10000 12500 15000];
% for testno = 1:10



for enabledchannel = 1:1

symbol          = [-1 -1 1 1 -1 -1 1 1  ...
                   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
                   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
                   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
                   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; % the symbol form
intervals       = 5000; % the number of intervals for which symbol+dumping will be repeated 
samplingRate    = 20e6; % the sampling rate in Hz for replay/recording

enablesave = true ;
enabledchannel;

%%




% ----- END OF USER INPUT ------
% DO NOT MODIFY CODE BELOW
close all
show_plots      = false;

% ----- do not alter this, for the time being ----
dumpPoints      = 0; % the number of dumping points after symbol
channels = 16;
plotOutSignal = false; % false if you do not want the output signal to be plotted
timeout_ms = 5000; % the replay timeout in ms for card 0
triggerOffset = 30; 

% ----- calculate the output signal train -----
[outDigData, outSignals, outSignal, outGate, outPol] = ...
    generateSignal(intervals, symbol, dumpPoints, channels, plotOutSignal);

outGateComplement = -1*repmat(outGate, 2*channels, 1) + 1;


%%
% %disable all channels above enabled channel 
% 
% if enabledchannel ~= 16
%     outDigData((2*enabledchannel)+1:32,:) = 1;
%     outSignals((enabledchannel+1):16,:) = 0;
% end

%%
%disable certain channels 
% 
% if enabledchannel ~= 16
%     outDigData((2*enabledchannel)+1:32,:) = 1;
%     outSignals((enabledchannel+1):16,:) = 0;
% end
% 
% if enabledchannel ~= 1
%     outDigData(1:(2*enabledchannel)-2,:) = 1;
%     outSignals(1:enabledchannel-1,:) = 0;
% end
%%



% ----- digital data -----
% calculate the muxed signal
outRAWData = bi2de(reshape(outDigData, 16, [])')';

starhubFound = false;

%  ----- assume 2 cards being present -----
for cardCount = 0 : 2

    % ----- init card and store info in cardInfos struct -----
    [success, cardInfo] = spcMInitCardByIdx (cardCount);

    if success == false
        % ----- check if at least one card is present -----
        if cardCount == 0
            spcMErrorMessageStdOut (cardInfo, 'Error: Could not open card\n', true);
            return;
        end

        break;
    else
        if bitand (cardInfo.featureMap, 32) | bitand (cardInfo.featureMap, 64)  % SPCM_FEAT_STARHUB5 = 32, SPCM_FEAT_STARHUB16 = 64
            starhubFound = true;
            fprintf ('Starhub found\n');
        end

        cardInfos (cardCount+1) = cardInfo;
    end
end        

% ----- printf info of all cards -----
for idx = 1 : cardCount 
    cardInfoText = spcMPrintCardInfo (cardInfos(idx));
    fprintf ('\nCard %d:\n', idx-1);
    fprintf (cardInfoText);
end

% ----- check if starhub is present-----
if starhubFound == false
    fprintf ('\nThere is no starhub in the system, this script can not run\n');
    starhubOk = false;
else
    starhubOk = true;  
end

% ----- the star hub is accessed under it's own handle -----
if starhubOk == true

    hSync = spcm_hOpen ('sync0');
    if hSync == 0
        fprintf ('\nCan not open starhub handle\n');
        starhubOk = false;
    end

    if starhubOk == true
        [errorCode, syncCards] = spcm_dwGetParam_i32 (hSync, 48990);  % 48990 = SPC_SYNC_READ_SYNCCOUNT

        % ----- show cable connection info -----
        fprintf ('\nStar-hub information:\n');
        fprintf ('Star-hub is connected with %d cards\n', syncCards);
        for idx = 0 : syncCards - 1
            [errorCode, cable] = spcm_dwGetParam_i32 (hSync, 49100 + idx);  % 49100 = SPC_SYNC_READ_CABLECON0
            fprintf ('  Card Idx %d (sn %05d) is', idx, cardInfos(idx+1).serialNumber);
            if cable ~= -1
                fprintf (' connected on cable %d\n', cable);
            else
                fprintf (' not connected with the star-hub\n');
            end
        end
        fprintf ('\n');

        % ----- first card will be set for replay -----
        % ----- set channel mask for max channels -----
        chMaskH = 0;    % this is always 0 for cards with 32 channels or less
        chMaskL = hex2dec ('FFFFFFFF'); % all bits to 1 -> all 32 channels enabled

        % ----- standard single, all channels, memsize=4*8*length(outSignal) -----    
        [success, cardInfos(1)] = spcMSetupModeRepStdSingle (cardInfos(1), chMaskH, chMaskL, 2*size(outDigData, 2));
%         [success, cardInfos(1)] = spcMSetupModeRepStdSingle (cardInfos(1), chMaskH, chMaskL,      67500);
        
        
        % ----- we try to set the samplerate on internal PLL, no clock output -----
        [success, cardInfos(1)] = spcMSetupClockPLL (cardInfos(1), samplingRate, 0);  % clock output : enable = 1, disable = 0

        fprintf ('  Replay Card %d: Sampling rate set to %.1f MHz\n', 0, cardInfos(1).setSamplerate / 1000000);

        % ----- set all output channel groups, no 110 ohm termination ----- 
        for i=0 : cardInfos(1).DIO.groups-1
            [success, cardInfos(1)] = spcMSetupDigitalOutput (cardInfos(1), i, 2, 0, 3300, 0);  % SPCM_STOPLVL_LOW = 2
        end

        % dwBytesPerSample: 2 = 16 bits, 1 = 8 bits 
        errorCode = spcm_dwSetRawData (cardInfos(1).hDrv, 0, length (outRAWData), outRAWData, 2);

        if (errorCode ~= 0)
            [success, cardInfos(1)] = spcMCheckSetError (errorCode, cardInfos(1));
            spcMErrorMessageStdOut (cardInfos(1), 'Error: spcm_dwSetData:\n\t', true);
            return;
        end

%         % ----- we'll start and wait until the card has finished or until a timeout occurs -----
%         errorCode = spcm_dwSetParam_i32 (cardInfos(1).hDrv, 295130, timeout_ms);  % 295130 = SPC_TIMEOUT
%         if (errorCode ~= 0)
%             [success, cardInfos(1)] = spcMCheckSetError (errorCode, cardInfos(1));
%             spcMErrorMessageStdOut (cardInfos(1), 'Error: spcm_dwSetParam_i32:\n\t', true);
%             return;
%         end


        % ----- set the remaining cards for acquisition -----
        for idx = 2 : cardCount

            % ----- standard single, all channels, memsize=2*length(data), posttrigge=2*length(data) -> pretrigger=0  -----    
            [success, cardInfos(idx)] = spcMSetupModeRecStdSingle (cardInfos(idx), chMaskH, chMaskL, 2*(size(outDigData,2)+triggerOffset) , 2*(size(outDigData,2)+triggerOffset));%512 * 1024, 512 * 1024);

            % ----- we try to set the samplerate to the samplingRate on internal PLL, no clock output -----
            [success, cardInfos(idx)] = spcMSetupClockPLL (cardInfos(idx), samplingRate, 0);  % clock output : enable = 1, disable = 0

            fprintf ('  Record Card %d: Sampling rate set to %.1f MHz\n', idx-1, cardInfos(idx).setSamplerate / 1000000);

            % ----- type dependent card setup -----
            switch cardInfos(idx).cardFunction
                   % ----- digital acquisition card setup (3 = DigitalIn, 5 = DigitalIO) -----
                case { 5 }
                    % ----- set all input channel groups, no 110 ohm termination ----- 
                    for i=0 : cardInfos(idx).DIO.groups-1
                        [success, cardInfos(idx)] = spcMSetupDigitalInput (cardInfos(idx), i, 0);
                    end
            end
        end

        % ------ Triggerring Setup -----

        % ----- 1st card is used as trigger master (un-comment the second line to have external trigger on card 0 -----
        [success, cardInfos(1)] = spcMSetupTrigSoftware (cardInfos(1), 0);  % trigger output : enable = 1, disable = 0
        %[success, cardInfos(1)] = spcMSetupTrigExternal (cardInfos(1), 1, 0, 0, 1, 0);  % 1 = SPC_TM_POS

        error = 0;
        syncEnableMask = bitshift (1, cardCount) - 1;
        syncClkMask = bitshift (1, (cardCount-1));

        % ----- sync setup, all card activated, last card is clock master -----
        error = error + spcm_dwSetParam_i32 (hSync, 49200, syncEnableMask);  % 49200 = SPC_SYNC_ENABLEMASK
        error = error + spcm_dwSetParam_i32 (hSync, 49220,    syncClkMask);  % 49220 = SPC_SYNC_CLKMASK

        % ----- start the card and wait for ready with the set time out -----
        error = error + spcm_dwSetParam_i32 (hSync, 295130, timeout_ms);  % 295130 = SPC_TIMEOUT

        if error == 0

            % ----- set command flags -----
            commandMask = bitor (4, 8);                %   M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER
            commandMask = bitor (commandMask, 16384);  % | M2CMD_CARD_WAITREADY

            fprintf ('\n  .... Acquisition started for all cards\n');

            error = spcm_dwSetParam_i32 (hSync, 100, commandMask);  % 100 = SPC_M2CMD 

            if errorCode == 263  % 263 = ERR_TIMEOUT   
                fprintf ('  .............................. Timeout');
            else
                fprintf ('  .... Acquisition successfully finished\n');

                analogDataIdx = 0;

                RAWData = zeros(1, cardInfos(2).setMemsize);

                for idx = 2 : cardCount

                    % ***** get digital input data *****
                    if cardInfos(idx).cardFunction == 5

                        % ----- get whole digital data in one multiplexed data block -----
                        % dwBytesPerSample: 2 = 16 bits (this cannot be changed without recoding spcMDemuxDigitalData)
                        [errorCode, RAWData] = spcm_dwGetRawData (cardInfos(idx).hDrv, 0, cardInfos(idx).setMemsize, 2);

                        % ----- demultiplex digital data (DigData (channelIndex, value)) -----
                        DigData = spcMDemuxDigitalData (RAWData, cardInfos(idx).setMemsize, cardInfos(idx).setChannels);
                        
                        % ----- convert the A/B bistream to analog signals 
                        AnalogSignals = ab_demux_signals(DigData);
                        
                        
                        % ----- cross-correlation                        
                        outSignalsNoDump = outSignals;
                        outSignalsNoDump( outSignalsNoDump == 2 ) = 0;
                        
                        DigDataExpanded = 2*DigData - int8(ones(size(DigData)));
                        DigDataGateComp = double(DigDataExpanded).*[ zeros(2*channels, triggerOffset) outGateComplement];
                        
                        mean1 = repmat(mean(outSignalsNoDump')', 1, size(outSignalsNoDump,2));
                        mean2 = repmat(mean(DigDataGateComp(1:channels, triggerOffset+1:end)')', 1, size(DigDataGateComp(1:channels, triggerOffset+1:end),2));
                        signalCorr = (outSignalsNoDump - mean1).*(DigDataGateComp(1:channels, triggerOffset+1:end) - mean2);
                        clear mean1 mean2
                        
                        meanOutSignals = repmat(mean(outSignalsNoDump')', 1, size(outSignalsNoDump,2));
                        meanAnalogSignals = repmat(mean(AnalogSignals(:, triggerOffset+1:end)')', 1, size(AnalogSignals(:, triggerOffset+1:end),2));
                        signalCorr = sum( (outSignalsNoDump - meanOutSignals).*(AnalogSignals(:, triggerOffset+1:end) - meanAnalogSignals), 2);
                        clear meanOutSignals meanAnalogSignals
                        
                        if (errorCode ~= 0)
                            [success, cardInfo(idx)] = spcMCheckSetError (errorCode, cardInfo(idx));
                            spcMErrorMessageStdOut (cardInfo(idx), 'Error: spcm_dwGetData:\n\t', true);
                            return;
                        end
                    end    
                end
            end
        end
    end
end

% ***** plot the data of channel 0 for each analog sync card *****
if show_plots
for idx = 2 : cardCount
    % ***** plot digital data (muxed) *****
%     figure;
%     spcMPlotDigitalData (DigData, cardInfos(idx).setChannels, size(DigData, 2));
%     xlim([triggerOffset size(DigData,2)-1]);
%     xlabel('sample');
%     ylabel('channel');
%     title('acquired digital signal');
    
    % ***** plot analog acquired data (demuxed) *****
%     figure;
%     plot(0:size(AnalogSignals,2) - 1, AnalogSignals + ...
%         reshape(repmat(6*[1:cardInfos(idx).setChannels/2],  size(AnalogSignals,2), 1),  size(AnalogSignals,2), cardInfos(idx).setChannels/2)');
%     xlim([triggerOffset size(AnalogSignals,2)-1]);
%     xlabel('sample');
%     ylabel('channel');
%     title('acquired analog signals');
    
    % ***** plot dump-less out signal and gating signals
    figure;
    plot(0:size(outSignalsNoDump,2) - 1, outSignalsNoDump + ...
        reshape(repmat(6*[1:cardInfos(idx).setChannels/2],  size(outSignalsNoDump,2), 1),  size(outSignalsNoDump,2), cardInfos(idx).setChannels/2)');
    hold on;
    plot(0:length(outGate) - 1, repmat(outGate, cardInfos(idx).setChannels/2, 1)  + ...
        reshape(repmat(6*[1:cardInfos(idx).setChannels/2],  length(outGate), 1),  length(outGate), cardInfos(idx).setChannels/2)', 'k');    
    xlim([triggerOffset size(AnalogSignals,2)-1]);
    xlabel('sample');
    ylabel('channel');
    title('Analog transmitted * gate');              %'acquired analog signals (gated, no dumping)');

    % ***** plot dump-less out signal and gating signals
    figure;
    plot(0:size(DigDataGateComp,2) - 1, DigDataGateComp + ...
        reshape(repmat(6*[1:cardInfos(idx).setChannels],  size(DigDataGateComp,2), 1),  size(DigDataGateComp,2), cardInfos(idx).setChannels)');
    hold on;
    temp1 = [zeros(2*channels, triggerOffset) outGateComplement];
    plot(0:size(temp1,2) - 1, temp1  + ...
        reshape(repmat(6*[1:cardInfos(idx).setChannels],  size(temp1,2), 1),  size(temp1,2), cardInfos(idx).setChannels)', 'k');    
    xlim([triggerOffset size(AnalogSignals,2)-1]);
    xlabel('sample');
    ylabel('channel');
    title('Received data +/-1 * gate_complement');          %'complemented analog signals (gated, no dumping)');

    
end
hold off;
end

% ***** close driver *****
if hSync ~= 0
    spcm_vClose (hSync);
end

for idx=1 : cardCount 
    spcMCloseCard (cardInfos(idx));
end
%%



%%%% START CODE HERE%%%%%%%%%%%%%%%%%%%%%%%%%

disp('  .... Assembling fmc')

x_max = 2000;
fmc = zeros(16,16,x_max);
for i = 1:16
    for j = 1:16
        switch j
            case [2 6 10 14 3 7 11 15]
                x = j+1;
            case [4 8 12 16]
                x = j-2;
            otherwise
                x = j;
        end
%         temp = xcorr(double(DigDataExpanded(x,:)), double(outSignalsNoDump(i,:)),x_max); 
%         temp = xcorr(double(DigDataExpanded(x,:)), double(my_sig(i,:)),x_max);  
%         temp = xcorr(double(DigDataExpanded(x,:)), double(signals(i,:)),x_max);  
        temp = xcorr(double(DigDataGateComp(x,:)), double(outSignalsNoDump(i,:)),x_max);  
        
        fmc(i,j,:)  = temp (x_max+2:end); 
    end
end

figure

subplot(2,1,1)
plot( (0:x_max-1) /samplingRate *1e6 , squeeze(fmc(1,1,:)))
% ylim([-.5 .5]*intervals/1000)
subplot(2,1,2)
plot( (0:x_max-1) /samplingRate *1e6 , squeeze(fmc(8,8,:)))
% ylim([-.5 .5]*intervals/1000)
xlabel('us')

% for test=1:16
%    figure
%    plot( (0:x_max-1) /samplingRate *1e6 , squeeze(fmc(test,test,:))) 
% end


 pathname = 'C:\DanielD\TestData\ProbeTest\BMode\' ;

 
if enablesave == true;
    display('  .... Saving Data')
    save(strcat(pathname,'bmodestep8'))
end
end
% end

