%% Wavelet-Fourier-Filtering for vertical stripes
% [Münch et al., 2012, Stripe and ring artifact removal with combined 
% wavelet — Fourier filtering]
% It utilizes a Fourier filtering approach in the following steps:
% 1. Decomposition: Wavelet decomposition
% 2. Fourier Filtering
%-------------------------------------------------------------------------
% Author: Niklas Rottmayer
% Date: 31.05.2023
%-------------------------------------------------------------------------
% Input:
% ima       -   corrupted image
% decnum    -   decomposition depths (dec > 0)
% wname     -   wavelet type (string)
% sigma     -   Gaussian damping parameter
%
% Output:
% nima      -   Result of stripe removal
%-------------------------------------------------------------------------


function [nima]=WFF(ima,decNum,wname,sigma)

    % Step 1: Wavelet decomposition
    for ii=1:decNum
        [ima,Ch{ii},Cv{ii},Cd{ii}]=dwt2(ima,wname);
    end

    % Step 2: Fourier filtering
    for ii=1:decNum
        % Fourier transform
        fCv=fftshift(fft(Cv{ii}));

        % Damping
        [my,mx]=size(fCv);
        damp=1-exp(-(-floor(my/2):-floor(my/2)+my-1).^2/(2*sigma^2));
        fCv=fCv.*repmat(damp',1,mx);

        % Inverse Fourier transform
        Cv{ii}=ifft(ifftshift(fCv));
    end

    % Step 3: Wavelet reconstruction
    nima=ima;
    for ii=decNum:-1:1
        nima=nima(1:size(Ch{ii},1),1:size(Ch{ii},2));
        nima=idwt2(nima,Ch{ii},Cv{ii},Cd{ii},wname);
    end
return