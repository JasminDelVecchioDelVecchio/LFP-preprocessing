function [cleaned_in] = ECG_ArtifactRemoval(in, sr, condition, eventQRS,time) 

% close all

% time = linspace(0, length(in)/sr , length(in));

for ii = 1:length(eventQRS)%making original signal templates based on qrs templates for decomposition
    start(ii) = round((time(eventQRS(ii))-0.20)*sr);
    stop(ii) = round((time(eventQRS(ii))+0.44)*sr);
end
rem1 = find(start<0);
rem2 = find(stop>length(in));
rem = [rem1,rem2];
start(:,rem) = [];
stop(:,rem) = [];

if length(start)>length(stop)
    start(end) = [];
elseif length(stop)>length(start)
    stop(end) = [];
end

for ich = 1:size(in,1)
    for ii = 1:length(start)
        seg = start(ii):stop(ii);
        if ii==1
            updated_template(ii,:) = in(ich,seg);
        else
            rem = length(seg) - length(start(1):stop(1));
            if rem > 0
                seg = start(ii):stop(ii)-rem;
                updated_template(ii,:) = in(ich,seg);
            elseif rem <0
                seg = start(ii):stop(ii) + abs(rem);
                updated_template(ii,:) = in(ich,seg);
            elseif rem == 0
                updated_template(ii,:) = in(ich,seg);
            end
        end
    end

    [U,S,V] = svd(updated_template', 0);%decomposing templates samples*time
    num_comp = 10;
    component = U(:,1:num_comp)*S(1:num_comp,1:num_comp);%plotting the first num_comp components
    componentOffset = .1;
    colors = [0.9,0,0.1; 0.8,0.7,1; 0.8,1,0.7; 0.8,1,0.2; 1,0.6,0;...
        1,0.8,0; 1,0.7,0.8; 0.3,0,0.9; 0.3,0.8,0.1; 0.6,0,0.1];
    figure;
    colororder(colors)
    % plot(component,'LineWidth',2.5); 
    hold on; % Ensure all plots are on the same axes

    for i = 1:size(component, 2)
        plot(component(:, i) + i * componentOffset, 'LineWidth', 2); % Add offset
    end
    
    xlabel('samples'); ylabel(' Voltage [mV]'); title('Component');
    legend('1','2','3','4','5','6','7','8','9','10','Location', 'best','Orientation','horizontal');
    grid on;
    figure;
    [pxx,f] = pwelch(in(ich,:),500,300,[],sr);%freq. spectrum of original signal
    plot(f,10*log10(pxx))

    while condition ==1
        if input('do you want to clean the signal?\n')%clean signal?
            bad_components = input('enter components to remove: \n');
            if bad_components~=0
                artefact = U(:,bad_components)*S(bad_components,:)*V';%calculate artefact based on components chosen
            end
            % note> this function is taken from brainstorm *it should be
            % open
            artefact_filtered = art_bandpass_hfilter(artefact',80,1,sr, 1, 0, 'filtfilt')';
            cleaned_template = updated_template - artefact_filtered';
            cleaned_signal(ich,:) = in(ich,:);
            for ii = 1:height(cleaned_template)
                cleaned_signal(ich,start(ii):start(ii)+length(cleaned_template(ii,:))-1) = cleaned_template(ii,:);%recreate the ecg removed signal
            end

        else
            disp('no components removed')
            cleaned_signal(ich,:) = in(ich,:);
        end
        %%
        time = (0:length(in)-1)./sr;
        figure;
        subplot(121)
        hold on
        plot(time, in(ich,:),'LineWidth',2,'Color','k') % plot original signal
        plot(time, cleaned_signal(ich,:),'LineWidth',2,'Color','b') % plot ecg removed signal
        xlabel('time (ms)');
        ylabel('Voltage (mV)');
        title('Uncleaned (+ECG) vs cleaned LFP')
        legend('uncleaned', 'cleaned')
        ax1 = gca; ax1.XGrid = 'on'; ax1.YGrid = 'on'; ax1.GridLineStyle = '--';  ax1.GridAlpha = 0.25;
        ax1.XMinorGrid = 'on'; ax1.YMinorGrid = 'on'; ax1.MinorGridLineStyle = '-.'; ax1.MinorGridAlpha = 0.15;
        %ax1.OuterPosition=[0 0 1 1];

        subplot(122)
        hold on
        [pxx,f] = pwelch(in(ich,:),sr,sr/2,[],sr);%freq. spectrum of original signal
        plot(f,10*log10(pxx),'LineWidth',2,'Color','k')
        % plot(f,pxx,'LineWidth',2,'Color','k')
        [pxx,f] = pwelch(cleaned_signal(ich,:),sr,sr/2,[],sr);%freq. spectrum of ecg removed signal
        plot(f,10*log10(pxx),'LineWidth',2,'Color','b')
        xlabel('Frequency [Hz]')
        ylabel('PSD [dB/Hz]')
        legend('uncleaned LFP(+ECG)', 'cleaned LFP')
        title('Frequency Spectrum')
        ax1 = gca; ax1.XGrid = 'on'; ax1.YGrid = 'on'; ax1.GridLineStyle = '--'; ax1.GridAlpha = 0.25;
        ax1.XMinorGrid = 'on'; ax1.YMinorGrid = 'on'; ax1.MinorGridLineStyle = '-.'; ax1.MinorGridAlpha = 0.15;


        satisfied = input('are you satisfied with the processing. 1 for yes, 0 for no.\n');%modify data structure?
        if satisfied == 1
            cleaned_in(ich,:) = cleaned_signal(ich,:);
            condition = 0;            
        end
    end
    condition = 1;
    close all;
end
end