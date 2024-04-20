%% Gabor Filter
% Generate the real part of Gabor filter. [Wikipedia - Gabor filter]
%-------------------------------------------------------------------------
% Author: Niklas Rottmayer
% Date: 02.06.2023
%-------------------------------------------------------------------------
% Input:
% sz        -   image size (array[2])
% lmb       -   wavelength of sinusoidal factor
% theta     -   orientation of the normal to the parallel stripes
% phi       -   phase offset
% sigma     -   standard deviation of Gaussian envelope
% gamma     -   spatial aspect ratio (ellipticity of support)
%
% Output:
% f         -   real part of Gabor filter
%-------------------------------------------------------------------------

function f = GaborFilter(sz,lmb,theta,phi,sigma,gamma)
    % centered in the image center
    [x,y] = meshgrid(-sz(2)/2+0.5:sz(2)/2-0.5,-sz(1)/2+0.5:sz(1)/2-0.5);
    f = exp(-( (x*cos(theta) + y*sin(theta)).^2 + gamma^2*(-x*sin(theta) + y*cos(theta)).^2 )./(2*sigma^2)) ...
            .* cos(2*pi*(x*cos(theta) + y*sin(theta))/lmb + phi);
end