function [lfp_embed, inds] = embedLFP(lfp, lags)
% input is a time x chan matrix

l = size(lfp, 1);

mi = min(lags);
ma = max(lags);

inds = (-mi+1:l-ma)';

lfp_embed = lfp(round(inds), :);

