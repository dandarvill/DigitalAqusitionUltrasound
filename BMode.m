clear all
close all

 %Path = 'G:\FYP\TestData\BMode\' ;
% file = '_repeat' ;
Path = 'H:\4th Year\';
%Path = 'C:\DanielD\TestData\ProbeTest\BMode\';
file = 'bmodestep';

%channel = [1 5 3 7 2 6 4 8 9 13 11 15 10 14 12 16]; %rearrange to match physical order
channel =  [1 2  4 5 6  7 8 ];

newsignal = zeros(8,5001);
shift = 1;
% for repeat = 1:3
    
    for j =6:6
        shift = 1;
        clearvars -except j Path file repeat shift newsignal channel
        load(strcat(Path, file,num2str(j)),'fmc', 'x_max', 'samplingRate') 
        
        figure
        for i = 1:7
            newsignal = squeeze(fmc(channel(i),channel(i),:));
            
            %Normalise and shift
            %newsignal =  newsignal./(max(newsignal*3));
            amp1(i) = max(newsignal(250:320));
            amp2(i) = max(newsignal(390:440));
            newsignal = newsignal + shift; 
                       
            plot( (0:x_max-1) /samplingRate *1e6, newsignal)                      
            hold on
            shift=shift+6000;
        end  
        
        xlabel('Time [\mus]')
        ylabel('Channel')
        xlim([10 40])
        ylim([0.1 7.9])
        title(strcat(num2str(j)))
        
        
        
    end   
% end

ratio = amp1./amp2; 
I = mat2gray(ratio);
figure
imshow(I,'InitialMagnification',10000)


%amp1 = amp1./max(amp1);
I = mat2gray(amp1);
figure
imshow(I,'InitialMagnification',10000)


% signal = signal/3;
% shift = 1;

% figure
% for i =1:8
% %     sd = 500;
% %     start = 50;
% 
%     %account for discrepencies in firing time
% %     while sd>70
% %         start=start+10;
% %         sd = std(signal(i,start:start+10));
% %     end
%     
% %     if start>2000
% %         start=1;
% %     end
%     
%     newsignal = squeeze(fmc(i,i,:));
%     
%     %Normalise
%     newsignal =  newsignal./(max(newsignal*3));
%    
%     newsignal = newsignal + shift;   
%     plot( (0:x_max-1) /samplingRate *1e6, newsignal)  
%     hold on
%     shift=shift+1;
%     
% end

% xlabel('Time [\mus]')
% ylabel('Channel')
% xlim([10 40])
% ylim([0.1 8.9])
% 
% title('B-Scan of High Side')

% xlim([0 70]) (0:x_max) /samplingRate *1e6,


%         [l, sz] = size(signal);
%         trig = 0;
                  
%         time = 1484e-3*((0:x_max) /samplingRate *1e6);
%         
%         plot(time, signal+shift)
%         hold on
%         title('B-Mode Imaging of constant slope')
%           
%         shift = shift+10000;

