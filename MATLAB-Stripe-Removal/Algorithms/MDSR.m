%% Multi-directional Stripe Remover
% [Liang et al., 2016, Stripe artifact elimination based on nonsubsampled
% contourlet transform for light sheet fluorescence microscopy]*
% * We restrict filtering to a subset of directional subimages by angle.
%
% It utilizes a Fourier filtering approach in the following steps:
% 1. Decomposition: Non-subsampled contourlet transform (NSCT)
% 2. Fourier Filtering
%-------------------------------------------------------------------------
% Author: Niklas Rottmayer
% Date: 31.05.2023
%-------------------------------------------------------------------------
% Input:
% img       -   corrupted image
% dir       -   number of decomposition directions (dir = 2^n)
% dec       -   decomposition depths (dec > 0)
% sigma     -   Gaussian damping parameters
% sigmaa    -   Gaussian adjustment parameter for sigma
% d         -   max. difference in directional angle for filtering
%
% Output:
% img       -   Result of stripe removal
%-------------------------------------------------------------------------
% Requirement:
% Nonsubsampled Contourlet Toolbox by Arthur Cunha (Version 1.0.0.0)
% Additionally, open the location of the add-on. Go into nsct_toolbox and
% type into the command line: 
%               mex atrousc.c
%               mex zconv2.c
%               mex zconv2S.c
% Now everything is installed correctly.

function img = MDSR(img,dir,dec,sigma,sigmaa,d)
    assert(d > 0,'d must be positive.')
    assert(sigma > 0,'sigma must be positive.')
    assert(dec >= 0,'The decomposition depths must be non-negative.')
    assert(dir > 0,'The number of decomposition directions must be positive.')
    % Angle difference
    angles = [pi/4 + pi*(1:2:dir)/(2*dir), 5*pi/4 - pi*(1:2:dir)/(2*dir)];
    dangles = min(abs(angles - pi/2),abs(angles - 3/2*pi));

    % Step 1: Apply nsct-decomposition (with default filters)
    C = nsctdec( img, repmat(round(log2(dir)),dec),'dmaxflat7', 'maxflat');

    % Step 2: Apply Fourier filtering to all vertical components
    % First entry <-> level (application to all "high-pass" levels, >1)
    % Second entry <-> direction (2 and 3 are vertical)
    for i = 1:dec
        % Damping
        for j = 1:dir
            sigmabar = exp(-dangles(j)^2/(2*sigmaa^2))*sigma;
            
            if dangles(j) <= d  
                % Fourier transform (fft + fftshift)
                fCv = fftshift(fft(C{i+1}{j}));
    
                % Damping
                [ny,nx]=size(fCv);
                damp=1-exp(-(-floor(ny/2):-floor(ny/2)+ny-1).^2/(2*sigmabar^2));
                fCv=fCv.*repmat(damp',1,nx);
    
                % Inverse Fourier transform (ifftshift + ifft)
                C{i+1}{j}=ifft(ifftshift(fCv));
            end
        end
    end

    % NSCT reconstruction (with default filters
    img = nsctrec( C, 'dmaxflat7', 'maxflat');
end