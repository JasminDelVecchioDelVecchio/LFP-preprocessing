function [dataFiltered] = FilterData(data, cutoff, fs, WL)
% This function filters marker tracks based on a lowpass filter and smoothing techniques described in Vissani et al. (2021)
% data - 2D matrix containing marker tracking data with dimensions (#samples x #coordinates).
% cutoff - cutoff frequency for the lowpass filter [Hz]
% fs - sampling frequency [Hz]
% WL - window length for smoothing using the Savitzky-Golay filter
   
    if rem(WL,2) == 0
        WL = WL-1;
    end
    
    flag = 0;
    
    % Most likely the number of coordinates is lower than the samples
    if size(data,1) < size(data,2)
        data = data';
        flag = 1;
    end

    % Interpolate data
    for y=1:size(data,2)
        indx_1=find(~isnan(data(:,y)));
        inan=find(isnan(data(:,y)));
        indx=find(~isnan(data(:,y)));
        tracksToInterpolate=data(indx,y);
        X=1:size(data,1);
        data(:,y)=interp1(X(indx),tracksToInterpolate',X,'spline');
        data(1:indx_1(1)-1,y)=NaN;
        data(indx_1(end)+1:end,y)=NaN;
    end
    
    for ii = 1:size(data,2)
        data2filt = [zeros(WL,1)+data(1,ii); data(:,ii); zeros(WL,1)+data(end,ii)];
        dataFiltLP = lowpass(data2filt, cutoff, fs);
        dataFiltSM = smooth(dataFiltLP, WL,  'sgolay');
        dataFiltered(:,ii) = dataFiltSM(WL+1:end-WL);
    end
    
    if flag 
        dataFiltered = dataFiltered';
    end
end

