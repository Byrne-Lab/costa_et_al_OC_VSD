function [filtered] = vsd_ellip(data)
% Simple Elliptic filter for VSD data. Fixed parameters.
[b,a] = ellip(2,0.1,40,[10 100]*2/1000);
filtered = filtfilt(b,a,data);
end

