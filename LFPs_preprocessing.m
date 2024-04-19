% Script for preprocessing GPi-LFP data in dystonia patients with DBS
clc
clear
close all

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
end

% Load subject data using 'bst_get' functions:
%   - 'bst_get(''Subject'', subj_name)' retrieves subject information
%   - 'bst_get(''StudyWithSubject'', sub.FileName)' retrieves associated studies
%   - 'listdlg' displays a dialog to select the appropriate study (resting state)

subj_name = 'sbj01'; % patient id

[sub,sub_idx] = bst_get('Subject', subj_name);

[studies,istudies] = bst_get('StudyWithSubject', sub.FileName);

[study_idx,~] = listdlg('ListString',{studies.Name}.');

resting_istudy = istudies(study_idx);
resting_data   = studies(study_idx).Data;
[~, Ntrial] = size(resting_data); % number of trials recorded on the same day
resting_ChannelMat = in_bst_channel(studies(study_idx).Channel.FileName);

% Loading options
LoadOptions.IgnoreBad      = 0;  % From raw files: ignore the bad segments
LoadOptions.ProcessName    = 'process_unknown';
LoadOptions.RemoveBaseline = 'no';
LoadOptions.UseSsp         = 1;

% for itrial = 1:Ntrial

itrial = 1;
    resting_DataMat    = in_bst_data(resting_data(itrial).FileName);

    [sMat, nSignals, iRows] = bst_process('LoadInputFile', resting_data(itrial).FileName, [], [], LoadOptions);
    [DataMat, matName] = in_bst(resting_data(itrial).FileName, [], 0);
    events=DataMat.F.events; % events labels

    Time = resting_DataMat.Time;
    Fs = 1/mean(diff(Time)); % sample frequency

    Data_struct = bst2eeglab(sMat.Data, sMat.Time, events, 1:size(sMat.Data,1), resting_ChannelMat, DataMat.ChannelFlag);

    % - Remove bad segments based on signal quality flags.
    [data,offset] = removeBAD_DB(Data_struct,1);

    % Get Channels labels
    ChannelLabels = {resting_ChannelMat.Channel.Name};

    % Get LFPs Channels
    GPI_Chs_idx = find(cellfun(@(x)~isempty(x),strfind({resting_ChannelMat.Channel.Type},'STN')));

    % Get Kinematic Channels
    KIN_Chs_idx = find(cellfun(@(x)~isempty(x),strfind({resting_ChannelMat.Channel.Type},'MISC')));

    % Get GPi-LFP time series
    LFP_signal_raw = data(GPI_Chs_idx,:);

    %% This section removes stimulation artifacts using the PARRM algorithm [Datin-van Rijn et al., 2021].
    % PARRM reconstructs and subtracts the periodic artifact from the LFP signal

    % Find the index of 'resting'
    FileName = resting_data.FileName;
    stim_on = strfind(FileName, 'on');

    % Display the result
    if stim_on

        stimRate = 180; % [Hz] Stimulation frequency (Starting point for period grid search) - it changes depending on the patient-specific optimal stimulation parameters
        Period = FindPeriodLFP(raw_LFP,[1,length(raw_LFP)-1],Fs/stimRate); % Find period

        % Define parameters for the PARRM filter:
        %   - perDist: Window size in the period domain for averaging samples.
        %   - winSize: Width of the window in samples for the PARRM filter.
        %   - skipSize: Number of samples to skip between windows.
        %   - winDir: 'both' uses samples from the past and future for filtering.

        perDist=0.01; % Period space window for averaging
        winSize=2000; % Width of the window in samples
        skipSize=20; % Samples to skip between windows
        winDir='both'; % Filter using past and future samples

        % Create the PARRM filter and apply it to remove the stimulation artifact
        PARRM = PeriodicFilter(Period,winSize,perDist,skipSize,winDir);

        LFP_signal = ((filter2(PARRM.',LFP_signal_raw','same')-LFP_signal_raw')./(1-filter2(PARRM.',ones(size(LFP_signal_raw')),'same')) + LFP_signal_raw')';
    else 
        LFP_signal = LFP_signal_raw;
    end

    %% This section removes cardiac artifacts from LFPs using Singular Value Decomposition (SVD).

    % - Get ECG time series
    ECG_Chs_idx = find(cellfun(@(x)~isempty(x),strfind({resting_ChannelMat.Channel.Name},'ES1')));
    ECG_signal = data(ECG_Chs_idx,:);

    % - Identify QRS peaks from the ECG channel (EMG probe placed on patient's chest)
    [eventQRS_smp] = detect_qrs(ECG_signal,1,length(ECG_signal),Fs); % identify qrs peaks

    % - Segment LFPs based on QRS peaks for artifact removal.
    % - Apply SVD to decompose LFP segments and remove ECG artifact components.
    time = (0:length(LFP_signal)-1)./Fs;
    [LFP_signal_noECG] = ECG_ArtifactRemoval(LFP_signal, Fs, 1, eventQRS_smp, time);

    %% This section processes kinematic data (if available) for head tremor analysis using head-mounted retroreflective markers tracked by an optoelectronic system.
    % Presence of data in `KIN_Chs_idx` indicates channels recording head kinematics.
    % This typically implies the patient exhibited dystonic head tremor during recording.

    if ~isempty(KIN_Chs_idx)

        % - Identify head markers (forehead, temples) based on channel labels (HEAD_MX: forehead; HEAD_RX: right temple, HEAD_LX: left temple).
        loi={'HEAD_MX_x_kin';'HEAD_RX_x_kin';'HEAD_LX_x_kin'};
        kin_labels = {Data_struct.chanlocs(KIN_Chs_idx).labels}';
        markers_tracks = Data_struct.data(KIN_Chs_idx,:)';

        [index_loi]=extract_markers(kin_labels',loi)';
        head=markers_tracks(:,index_loi(1):index_loi(1)+2);
        head_rx=markers_tracks(:,index_loi(2):index_loi(2)+2);
        head_lx=markers_tracks(:,index_loi(3):index_loi(3)+2);

        % - Apply low-pass filtering and smoothing of the estimated head
        % angular velocity (Vissani et al., 2021)
        head = FilterData(head, 8, Fs, 30);
        head_rx = FilterData(head_rx, 8, Fs, 30);
        head_lx = FilterData(head_lx, 8, Fs, 30);


        % - Calculate head tilt, obliquity, and rotation angles. Head
        % reference point is defined as the midpoint (med) between temples
        % (assuming head symmetry)
        med=(head_rx+head_lx)/2;
        tilt=[(head(:,2)-med(:,2)),(med(:,1)-head(:,1))*-1];
        obl=[(head_rx(:,2)-head_lx(:,2)),(head_rx(:,3)-head_lx(:,3))*-1];
        rot=(head_rx(:,[1,3])-head_lx(:,[1,3]));

        tilt_angle=(atand(tilt(:,1)./tilt(:,2)));
        rotation_angle=(atand(rot(:,1)./rot(:,2)));
        obliquity_angle=(atand(obl(:,1)./obl(:,2)));

        head_angular_positions=[tilt_angle rotation_angle obliquity_angle];

        % - Compute head angular velocity using MarkersVelocity function for tremor assessment.
        h = 1; % incremental step
        vel = MarkersVelocity(head_angular_positions, h, Fs);

        % % - Apply high-pass filter to the estimated head angular velocity
        ord = 2; fc = 1; fn = Fs/2;
        [b,a] = butter(ord,fc/fn,'high');
        vel_filt = filtfilt(b,a,vel);
        vel = vel_filt;

        % - Remove bad segments based on signal quality flags.
        bad_segments = Data_struct.bad_segment;
        vel(bad_segments,:) = [];

    end

    % - Regression analysis to remove head tremor artifact from LFPs (Del Vecchio Del Vecchio et al., 2023)
    lags = linspace(-Fs,Fs,51); %2-sec window with K = 51 equally spaced shifts for temporal embedding
    % - Embed head angular velocity (vel) with time lags for regression analysis
    vel_embedded = embedHeadAngularVelocity(vel, lags)';
    % - Embed LFP signal for regression
    lfp_embedded = embedLFP(LFP_signal_noECG', lags)';


    % - Ensure equal time window for regression
    if size(vel_embedded,2)>size(lfp_embedded,2)
        % - Adjust shorter dimension to match the longer one (trimming the last element)
        vel_embedded(:,end)=[];
    elseif size(lfp_embedded,2)>size(vel_embedded,2)
        lfp_embedded(:,end)=[];
    end

    % - Regress out head tremor influence (embedded vel) from LFPs
    LFP_signal_cleaned = padHeadAngularVelocity((lfp_embedded - (lfp_embedded/vel_embedded)*vel_embedded)', lags, 0)';

    % - Compare PSDs of LFP signals after cleaning steps (Welch's method)
    % Compute PSDs of LFP with ECG and head tremor artifact removal
    [PSD_LFP_raw, f] = pwelch(LFP_signal', Fs, Fs/2, [], Fs);
    [PSD_LFP_noECG, ~] = pwelch(LFP_signal_noECG', Fs, Fs/2, [], Fs);
    [PSD_LFP_cleaned, ~] = pwelch(LFP_signal_cleaned', Fs, Fs/2, [], Fs);

    % Plot the comparison (optional)
    % figure;
    % plot(f, PSD_LFP_raw, 'LineWidth',1); hold on;
    % plot(f, PSD_LFP_noECG, 'LineWidth',1);
    % plot(f, PSD_LFP_cleaned, 'LineWidth',1);
    % legend('Raw', 'After ECG Removal', 'Final Cleaned');
    % xlabel('Frequency (Hz)'); ylabel('Power Spectral Density');
    % title('Comparison of LFP PSDs');

    % - Save cleaned LFP data for further analysis
    clearvars -except subj_name resting_data resting_ChannelMat itrial LoadOptions
% end
