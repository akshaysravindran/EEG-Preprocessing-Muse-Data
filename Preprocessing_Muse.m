%% MUSE EEG Pre-Processing Tutorial
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                        %
%                  EEG Preprocessing Tutorial Code                       %
%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Summary: From zero to hero in EEG pre-processing;  
% The code will help you navigate how to clean an EEG data collected using
% Muse system. 
% The code will focus on the following sections
%
% 0)  Load the required libraries
% i)  Load the EEG,chanlocs into EEGlab format
% ii) Downsample if needed
% iii)High Pass filter (0.1-2 Hz: Depending on noise level & study)
% iv)  Low pass filter (depending on the study)
% v) Remove Bad segments (Visually inspecting)
% vi)Remove Bad channels (Visual Inspection with Muse) 
% Vii) Remove burst artifacts using ASR 

% Additional Commands
% xi) Notch filter  
%
%
% Here by default EEG.data has the dimension channels x samples; If using
% custom data/matlab variable, ensure the dimension is correct
%
% Author: Akshay Sujatha Ravindran
% email: akshay dot s dot ravindran at gmail dot com
% Jan 20th 2021

plot_check = 1;
%% ia. Add EEGlab to the path 
% addpath('\eeglab2019');% Replace this line with the respective folder in your directory
% % To add permanently, go to home> set path > add folder > <Add the folder corresponding to EEGlab>
% eeglab % run this to load all dependencies and subfolders

% **************************************************************************************************************************************************
%% ib. Load the EEG variable
% Load the EEG data in .mat format

name        ='Data' % Name of the file
load([name,'.mat']);
% Assumes data is converted using museplayer. If not replace the next lines
EEGraw      = IXDATA.raw.eeg.data'; % EEG: length of data x number of channels
Acc_raw     = IXDATA.raw.acc.data'; % Accelerometer: length of data x number of channels

% Add the file into EEGLab format
EEG         = pop_importdata('dataformat','array','nbchan',size(EEGraw,1),'data','EEGraw','srate',Fs,'pnts',length(EEGraw));
%**************************************************************************************************************************************************

%% ic. Load the channels
% Add EEG channel location to the EEG structure variable
% Template: EEG_output_variable.chanlocs = pop_chanedit(EEG_output_variable.chanlocs, 'load',{ 'directory\EEG_channel_location_file.ced', 'filetype', 'autodetect'});
EEG.chanlocs=pop_chanedit(EEG.chanlocs, 'load',{ 'eggs.ced', 'filetype', 'autodetect'});

% Checks: Add the directory name if the channel location file is stored in different folder
% **************************************************************************************************************************************************
%% 2 Downsample if required
% Downsample into a smaller sampling frequency

% Double click the EEG structure variable in the workspace and check the srate variable inside the EEG structure; 
% Otherwise type the following in command window
% >> disp(EEG.srate);

% Template: EEG_output_variable = pop_resample(EEG_input_variable,Required_sampling_frequency); 

EEG         = pop_resample(EEG,110); % resample from 220 Hz to 110 Hz 
rawEEG      = EEG; % For plotting purpose

% >> Check <<
% Make sure the resampled frequency is atleast twice the frequency
% of interest Nyquist criterion; If looking at 50 Hz component, sampling
% rate should be 100+ Hz

% **************************************************************************************************************************************************
%% 3 High Pass filter
% Remove slow drift from the EEG signal


% Template: [filter coefficients] = butter(filter_order ,Sampling_rate/2, 'Type_of_filter'); 

Fc_hpf      = 1; % Cut off frequency
filt_order  = 2;     
[a,b]       = butter(filt_order,Fc_hpf/(EEG.srate/2),'high');  % Create the butterworth filter coefficients
EEG.data    = filtfilt(a,b,double(EEG.data)')'; 

hpfEEG      = EEG; % For plotting purpose
disp('High pass filtered')


% **************************************************************************************************************************************************
%% 4 Low Pass filter 
% Can replace Step 3 & 4 with band pass instead of HPF and LPF being applied separtely;\
% Using separate if one do not wish to use LPF, just dont run this cell.

% Remove higher frequency contents from the data; Minimize muscular
% artifacts

% Template: [filter coefficients] = butter(filter_order ,Sampling_rate/2, 'Type_of_filter'); 

Fc_lpf      = 15; % Cut off frequency
[a,b]       = butter(filt_order,Fc_lpf/(EEG.srate/2),'low');  % Create the butterworth filter coefficients
EEG.data    = filtfilt(a,b,double(EEG.data)')';

lpfEEG      = EEG; % For plotting purpose
disp('Low pass filtered')
%% %% 5 Remove Segments  (If needed)
% Redraw the eeglab GUI and remove segments

eeglab redraw
% Can remove through scripts; Preferable to check visually and remove using the GUI: Tutorial inside
% the slides

% Remove EEG samples between the points inside the bracket; use start and
% end sample number; If multiple segments are to be removed, use semi colon
% and add the next entry
% Template: Remove two windows corresponding to two pairs of begin and end
% EEG_output_variable = eeg_eegrej(EEG_input_variable, [begin_sample_1  end_sample_1; begin_sample_2  end_sample_2]);

if false % this is just to show the format for running the code in script, samples removed are random so dont execute this
    EEG     = eeg_eegrej( EEG_LPF, [1 , 1000; 150000 158000]); % segment between 1:1000 and 150000:158000 is removed   
end



% **************************************************************************************************************************************************

%% 6 Remove channels (If needed)

eegplot(EEG.data,'srate',EEG.srate,'winlength',20)
if false    
    EEG     = pop_select(EEG,'nochannel',[1]); % remove channels mentioned in the bracket
end

% **************************************************************************************************************************************************

%% Sanity Check Required
% >> Sanity Check 1 <<
% Plot the PSD always to check if the filter is doing what you expect

nw          = 4; % time-bandwidth product parameter
Ch_Name     = {EEG.chanlocs.labels}; % If no chanlocs {'TP9', {'AF7'} ,   {'AF8'} ,   {'TP10'}}
ch_to_plot  = 1; % Channel number to plot
st_sample   = EEG.pnts-100*EEG.srate; % Start sample for plotting
end_sample  = EEG.pnts-40*EEG.srate; % End sample for plotting
freq        = [0:0.1:EEG.srate/2]; % frequencies of interest   
pxx_original= pmtm(rawEEG.data(:,st_sample:end_sample)',nw,freq,EEG.srate);         % Compute PSD for raw EEG        
pxx_HPF     = pmtm(hpfEEG.data(:,st_sample:end_sample)',nw,freq,EEG.srate);         % Compute PSD for High pass filtered EEG  
pxx_LPF     = pmtm(lpfEEG.data(:,st_sample:end_sample)',nw,freq,EEG.srate);         % Compute PSD or Low pass filtered EEG    


% Plot all the PSDs together
figure
plot(freq, 10.*log10(pxx_original(:,ch_to_plot)),'linewidth',1.5)
hold on
plot(freq, 10.*log10(pxx_HPF(:,ch_to_plot)),'linewidth',1.5)
plot(freq, 10.*log10(pxx_LPF(:,ch_to_plot)),'linewidth',1.5) 
line([Fc_hpf Fc_hpf],[-60,45],'linestyle',':','linewidth',1.5,'color','k')
line([Fc_lpf Fc_lpf],[-60,45],'linestyle','--','linewidth',1.5,'color','k')

ylim([-60,45])
title([Ch_Name{ch_to_plot},' Channel PSD Comparison']) 
legend('Original','HPF','LPF','HP cutoff','LP cutoff')
ylabel('PSD (in dB)')
xlabel('Frequency (Hz)')
set(gca,'FontName','Times New Roman','fontsize',14)
set(gca, 'XColor', [0 0 0], 'YColor', [0 0 0])   
% **************************************************************************************************************************************************
EEG_orig=EEG;
%% 7. Artifact Subspace Reconstruction
% There are many additiona features inside, one imp being window size. Fine tuning
% those might provide better cleaning. Window size is currently set at half
% the sampling rate. ASR requires the window size to be
% big enough to capture the artifact, try changing it.

EEG         = clean_artifacts(EEG_orig, 'FlatlineCriterion','off',...  % Remove flatline channels
'Highpass','off',... % Keep it off as we already performed high pass filter
'ChannelCriterion','off',... % Keep it off if we are removing channels manually, Do not use with Muse
'LineNoiseCriterion', 'off',... % Keep it off if dont want to remove channels based on line noise criterion
'BurstCriterion',5, ... % Standard deviation cutoff for removal of bursts; Try modifying this and check if is not cleaning too much/too little
'WindowCriterion','off'); % Keep it off if you dont want to remove too noisy windows
% vis_artifacts(EEG_ASR,EEG_CR) : See all channels before and after ASR
asrEEG      = EEG;

%%
% %>> Sanity Check 2 <<
 ft_hastoolbox('brewermap', 1);
% Plot to check how the ASR is performing
N           = 1; % Start sample for plotting
N1          = EEG.pnts; % End sample for plotting
ch          = 1;
figure;
h(1)        = subplot(3,1,1)
plot([N:N1]./EEG.srate, EEG_orig.data(ch,N:N1));
hold on
plot([N:N1]./EEG.srate, asrEEG.data(ch,N:N1));
xlim([N N1]./EEG.srate)

xticklabels('')
ylabel('Mag (uV)' )
title('Channel TP9; ASR Threshold = 50 SD')
legend('Before ASR','After ASR','location','northwest')
legend('boxoff')
legend('location','southeast')
set(gca,'FontName','Times New Roman','fontsize',14)
set(gca, 'XColor', [0 0 0], 'YColor', [0 0 0])  



% Plot the spectrogram for LPF data
% Compute the spectrogram of the data 
[s,f,t,p]=spectrogram(lpfEEG.data(ch,N:N1),EEG.srate,EEG.srate-1,freq,EEG.srate);
h(2)        = subplot(3,1,2)
A=10.*(log10(p));
imagesc(t,(f),(A))
set(gca,'YDir','normal');

colormap(flipud(brewermap(100,'RdBu')));
caxis([-(max(A(:))),abs(max(A(:)))])
title('Before ASR')
ylabel('Freq(Hz)' )
set(gca,'FontName','Times New Roman','fontsize',14)
set(gca, 'XColor', [0 0 0], 'YColor', [0 0 0])  
xticklabels('')
L=colorbar('location','southoutside')
ylabel(L, 'dB')

% Plot the spectrogram for ASR cleaned data
h(3)        = subplot(3,1,3)
[s,f,t,p]   = spectrogram(asrEEG.data(ch,N:N1),EEG.srate,EEG.srate-1,freq,EEG.srate);
A1=10.*(log10(p));
imagesc(t,(f),(A1))
set(gca,'YDir','normal');
caxis([-(max(A(:))),abs(max(A(:)))])
ylabel('Freq(Hz)' )
xlabel('Time (s)')
title('After ASR')
linkaxes(h,'x')
colormap(flipud(brewermap(100,'RdBu')));
set(gca,'FontName','Times New Roman','fontsize',14)
set(gca, 'XColor', [0 0 0], 'YColor', [0 0 0])  
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 8, 8], 'PaperUnits', 'Inches', 'PaperSize', [8.5, 11])

% pop_saveset(EEG,'filename' ,'EEG_cleanVF','filepath','') % save final cleaned data
% **************************************************************************************************************************************************

%% 8 . Additional useful scripts
if false % change to true if you want to run this. Add the parts to the appropriate parts
    
    % 1) Line Noise Removal 
    % Remove 60 Hz line noise using notch filter
    Q_fact   = 20;
    Fc       = 60;
    wo       = Fc/(EEG.srate/2);  
    bw = wo/Q_fact;
    [b,a]    = iirnotch(wo,bw); 
    EEG.data = filtfilt(b,a,double(EEG.data'))'; 
    
    

end

        
