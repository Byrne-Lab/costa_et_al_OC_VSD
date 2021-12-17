function [pixels,kern_center,kernel_size,kernpos]=readdet_var(det)
% Function for reading kernel coordinates from det VARIABLE. Modified from
% readdet function, which reads kernel coordinates from .det file and
% generates the variable ("det") used as input here.
% "det" --  Nx1 double variable containing linear indexes of pixels
%           belonging to each kernel. Kernels are separated by zeros. This
%           variable is created by the "readdet" function.
 
% Parameters

pixSz = 128; % Size of the acquired image on which kernels were drawn. Currently assumes X and Y dimensions are the same.

pixels = mod(det,pixSz);
pixels(:,2) = ceil(det/pixSz);
pixels(pixels(:,1)==0,1) = pixSz;

kernpos = find(~det);
kernel_size = diff(kernpos)-1;
kernpos(end) = [];
kern_center = zeros(length(kernpos),2);
for a=1:length(kernpos)
    if a<length(kernpos) 
        B = pixels(kernpos(a)+1:kernpos(a+1)-1,:);
    else
        B = pixels(kernpos(a)+1:end-1,:);
    end
    kern_center(a,:) = mean(B);
end

end


