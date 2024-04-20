% RMS Contrast
% Root mean square (RMS) contrast does not depend on the angular frequency content or the spatial distribution of contrast 
% in the image. RMS contrast is defined as the standard deviation of the pixel intensities

function contrast = contrastMetric(im)
im2 = im2double(im);
N = size(im);
m= N(1,1); n=N(1,2);
I_hat = mean2(im2);
rms = 0;
for i=1:m-1
    for j=1: n-1
       rms = rms +(im2(i,j) - I_hat)^2; 
    end
end

contrast = sqrt((1/(m*n))*rms);

