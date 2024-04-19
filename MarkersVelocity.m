function [v] = markersVelocity(tracks, h,  frequency)
% This function calculates the velocity of markers in a tracking data structure using first-order forward differentiation.
% tracks - A 2D matrix containing tracking data with dimensions (#samples x #coordinates)
% h - incremental step size (integer)
% frequency - sampling frequency [Hz]
    
    flag = 0;
    
    % Most likely the number of coordinates is lower than the samples
    if size(tracks,1) < size(tracks,2)
        tracks = tracks';
        flag = 1;
    end
    
    tracks = [tracks; zeros(h,size(tracks,2))+tracks(end,:)];
    
    nFRM    = size(tracks,1);
    nTracks = size(tracks,2);
    
    v = zeros(nFRM-h, nTracks);
    
    for ii = 2:nFRM-h
        for jj = 1:nTracks
            v(ii,jj) = (tracks(ii+h,jj) - tracks(ii-h,jj))*0.5*frequency; %error of order 1/fs^2
        end
    end 
    
    if flag 
        v = v';
    end
end

