clc
clear
% close all

% Check if Brainstorm toolbox is accessible using 'bst_get' function
if exist('bst_get') == 2
    disp('Brainstorm is installed and available.');
else
    disp('Brainstorm is not installed or not available');
end

if exist('eeg_checkset') == 2
    disp('EEGLAB is installed and available.');
else
    disp('EEGLAB is not installed or not available');
    % C:\Users\jasmi\Lab Isaias Dropbox\Team\WUMIGE_shared\WUMIGE_shared_code\eeglab14_1_2b
end

% Subject
subj_name = 'PW_TB89';

% Right (rx) and Left (lx) suffixes (for kinematic events)
rx_lx_suffix = {'rx','lx'};

% Frequency vector [Hz]
freqs = 1:95;

% Load subject data using 'bst_get' functions:
%   - 'bst_get(''Subject'', subj_name)' retrieves subject information
%   - 'bst_get(''StudyWithSubject'', sub.FileName)' retrieves associated studies
%   - 'listdlg' displays a dialog to select the appropriate study (resting state)

fprintf(['Subject ' subj_name])

[sub, sub_idx] = bst_get('Subject', subj_name);

[studies,istudies] = bst_get('StudyWithSubject',   sub.FileName);

[study_idx,~] = listdlg('ListString',{studies.Name}.');

walking_data   = studies(study_idx).Data;

ChannelMat = in_bst_channel(studies(study_idx).Channel.FileName);

% Loading options
LoadOptions.IgnoreBad      = 0;  % From raw files: ignore the bad segments
LoadOptions.ProcessName    = 'process_unknown';
LoadOptions.RemoveBaseline = 'no';
LoadOptions.UseSsp         = 1;

for iFile = 1:length(walking_data) % all recorded trials
    % ===== GET CHANNEL FILE =====
    % Load channel file
    % Check for different channel structures in FilesA
    % Get channels to process

    [sMat, nSignals, iRows] = bst_process('LoadInputFile', walking_data(iFile).FileName, [], [], LoadOptions);
    [DataMat, matName] = in_bst(walking_data(iFile).FileName, [], 0);

    ProtocolInfo = bst_get('ProtocolInfo');

    MatFile = load(bst_fullfile(ProtocolInfo.STUDIES, walking_data(iFile).FileName));

    events = DataMat.F.events;
    %     clear DataMat;

    iChannels = channel_find(ChannelMat.Channel, {'STN'});
    excludeChannels = setdiff(1:size(sMat.Data,1),iChannels);

    EEG = bst2eeglab(sMat.Data,sMat.Time,events,1:size(sMat.Data,1),ChannelMat,DataMat.ChannelFlag);

    if iFile==1
        EEGmerged = pop_select(EEG, 'nochannel', excludeChannels);
        EEGmerged.badchan(excludeChannels) = [];
    else
        EEG = pop_select(EEG, 'nochannel', excludeChannels);
        EEG.badchan(excludeChannels) = [];
        EEGmerged = pop_mergeset(EEGmerged,EEG);
        EEGmerged.bad_segment = cat(2,EEGmerged.bad_segment,EEG.bad_segment);
        EEGmerged.badchan = EEGmerged.badchan | EEG.badchan;
    end

end

clear EEG;
str = 'HC_LX';

addpath 'C:\Users\jasmi\Lab Isaias Dropbox\Team\WUMIGE_shared\WUMIGE_shared_projects\Dystonia_jdvdv\EEG_analysis\CleaningPipeline\EEG-preprocessing\scripts'
addpath 'C:\Users\jasmi\Lab Isaias Dropbox\Team\WUMIGE_shared\WUMIGE_shared_code\eeglab14_1_2b\plugins\firfilt1.6.2'

[EEGmerged, com, b] = mypop_firws(EEGmerged, [], 80); %filtering data
EEGepoched_hc = pop_epoch(EEGmerged,{str},[-0.575 0.575]); %average stride length = 1.15 sec

% Visually check gait artifacts presence in all recorded LFP channels
fig1 = figure();
for jj=1:size(EEGepoched_hc.data,1) % number of LFP channels
    ax(jj) = subplot(2,2,jj);
    data_STN=reshape(EEGepoched_hc.data(jj,:,:),[size(EEGepoched_hc.data,2),size(EEGepoched_hc.data,3)])';
    imagesc(EEGepoched_hc.times,1:EEGepoched_hc.trials,data_STN);
    xline(0,'Color','w','Linewidth',2,'LineStyle','--');
    ind=regexp(EEGmerged.chanlocs(1,jj).labels,'seeg');
    title([strrep(EEGmerged.chanlocs(1,jj).labels(1:ind-2),'_',' '),' - ', strrep(str,'_',' ')]);
    xlabel('time (ms)')
    ylabel(['EPOCHS (centred at ', strrep(str,'_',' '),')'])
    indstr=max(regexp(EEGmerged.chanlocs(1,jj).labels(1:ind-2),'_'));
end
linkaxes(ax,'x')
sgtitle(subj_name{subjidx},'Interpreter','none')
clear ax

% average of each window of 1.15sec to find mean artifact template
artifact_template = mean(EEGepoched_hc.data,3);

% fitting of artifact template
artifact_template_fitted = [];

for istride = 1:size(EEGepoched_hc.data,3)
    for ich = 1:size(EEGepoched_hc.data,1)
        B = regress(EEGepoched_hc.data(ich,:,istride)',...
            [ones(size(artifact_template(ich,:)')), artifact_template(ich,:)']);
        artifact_template_fitted(:,:,ich,istride) = B(2)*artifact_template(ich,:)' + B(1);
    end
end

artifact_template_fitted = squeeze(artifact_template_fitted);
artifact_template_fitted = permute(artifact_template_fitted,[2,1,3]);
% removal of artifact template
Fs = 250;
fc = 1;
FWHM = 3;

data_cleaned = EEGepoched_hc.data - artifact_template_fitted;

% Data in Frequency Domain -------------------------------------- %
window   = Fs; % 250 samples
noverlap = window/2; % 50% overlap
nfft = window;


for istride = 1:size(data_cleaned,3)
    [pxx_cleaned(:,:,istride), fxx] = pwelch(data_cleaned(:,:,istride)', window, noverlap, nfft, Fs);
    [pxx(:,:,istride), ~] = pwelch(EEGepoched_hc.data(:,:,istride)', window, noverlap,nfft, Fs);
end


% Plot comparison between uncleaned and cleaned LFPs
% Gait artefact seems to be affecting mostly low-frequency oscillations 
figure();
plot(fxx, mean(pxx, 3, 'omitmissing'), 'LineWidth', 2); grid on; hold on;
plot(fxx, mean(pxx_cleaned, 3, 'omitmissing'), 'LineWidth', 2);
legend({'GPi-L raw', 'GPi-R raw', 'GPi-L cleaned', 'GPi-R cleaned'});
