%% General Stripe Remover (3D)
% This algorithms performs a variational stripe removal approach.
% It utilizes the following objective functional:
% argmin mu1*||D_(x,y,z) u||_(2,1) + ||D_y s||_1 + mu2*||s||_1 + iota_[0,1](u).
%  u+s=F
% Optimizer: PDHGMp [Burger et al., 2014, First Order Algorithms in
% Variational Image Processing]
%-------------------------------------------------------------------------
% Author: Niklas Rottmayer
% Date: 31.05.2023
%-------------------------------------------------------------------------
% Input: 
% F             -   corrupted image
% iterations    -   number of steps to optimize
% mu            -   weighting parameters (array[2])
% proj          -   project onto [0,1] (true/false)
% resz          -   ratio of resolutions z-axis to x and y (<=1) 
%                   D_z = resz*(F(:,:,2:end) - F)
% normalize     -   normalize the input image to [0,1] (true/false) 
% GPU           -   use GPU version (true/false)
% verbose       -   print process information (true/false)
%
% Output:
% u             -   Result of stripe removal
%-------------------------------------------------------------------------
% Comment: Suggested weighting parameters are:
%           - [10,0.1] generally good 
%           - [15,0.5] for more difficult stripes (e.g. slightly oblique)
%           - [12,0.2] for heavy stripes and [10,0.1] is not enough
%           - [ 7,0.1] in case [10,0.1] does too much (e.g. smoothing)

function u = GSR3D(F,iterations,mu,proj,resz,normalize,GPU,verbose)
    % Set default values
    if nargin < 8
        verbose = 1;
    end
    if nargin < 7
        GPU = 1;
    end
    if nargin < 6
        normalize = 1;
    end
    if nargin < 5
        resz = 1;
    end
    if nargin < 4
        proj = 1;
    end
    if nargin < 3
        mu = [10,0.1];
    end
    if nargin < 2
        iterations = 1000;
    end

    % Sanity checks:
    assert(length(size(F)) == 3,'The input image must be 3-dimensional.')
    assert(iterations >= 0,'The number of iterations must be non-negative.')
    iterations = ceil(iterations); % -> Safety conversion to integer!
    assert(length(mu) == 2,'mu must be an array of length 2.')
    assert(all(mu > 0),'mu must contain positive values.')
    assert(proj == 0 || proj == 1,'proj must be boolean (true/false or 1/0)')
    assert(resz > 0, 'resz must be positive.')
    if resz > 1
        warning('resz is larger than 1. This may cause problems.')
    end
    if verbose
        fprintf('Intializing 3D Stripe Removal\n... please wait ...')
    end
    % Initialization
    tau = 0.28; sigma = tau; % ~1/sqrt(12) (should be correct)
    lambda = 30; % Rescaling of optimization functional (no explanation)
    [nx,ny,nz] = size(F);
    if normalize
        immax = max(F,[],'all');
        immin = min(F,[],'all');
        F = (F - immin)/(immax - immin);
    end
    
    % Helper variables
    if GPU
        b1x = gpuArray(zeros([nx,ny,nz])); 
        F = gpuArray(F);
    else 
        b1x = zeros([nx,ny,nz]);
    end
    b1xbar = b1x; b1y = b1x; b1ybar = b1x; b1z = b1x; b1zbar = b1x; 
    b2 = b1x; b2bar = b1x; b3 = b1x; b3bar = b1x;
    
    % Step 0:
    u = F;
    s = b1x;
    
    % Main loop:
    k = 1;
    stopping = inf;

    while k <= iterations && stopping > 10^(-5)
        if verbose
            fprintf('\rIteration: %d / %d', k, iterations)
        end
        
        % Part 1: Update u and s
        % s1x = D_x^T b1xbar  
        s1x = b1xbar(:,[1 1:end-1],:) - b1xbar;
        s1x(:,1,:) = -b1xbar(:,1,:);
        s1x(:,end,:) = b1xbar(:,end-1,:);

        % s1y = D_y^T b1ybar
        s1y = b1ybar([1 1:end-1],:,:) - b1ybar;
        s1y(1,:,:) = -b1ybar(1,:,:);
        s1y(end,:,:) = b1ybar(end-1,:,:);

        % s1z = D_z^T b1zbar
        s1z = resz.*(b1zbar(:,:,[1 1:end-1]) - b1zbar);
        s1z(:,:,1) = -resz.*b1zbar(:,:,1);
        s1z(:,:,end) = resz.*b1zbar(:,:,end-1);

        % s2 = D_y^T b2bar
        s2 = b2bar([1 1:end-1],:,:) - b2bar;
        s2(1,:,:) = -b2bar(1,:,:);
        s2(end,:,:) = b2bar(end-1,:,:);

        u = u - tau * sigma * (s1x + s1z + s1y);
        s = s - tau * sigma * (s2 + b3bar);

        % Re-projection onto u+s=F (and u in [0,1])
        temp = F - s - u;
        u = u + 1/2 * temp;
        s = s + 1/2 * temp;
        if proj
            s = s + (u < 0).*u + (u > 1).*(u-1);
            u = max(min(u,1),0);
        end

        % Part 2: Updating helper variables (b and bbar)
        b1xbar = b1x; b1ybar = b1y; b1zbar = b1z; b2bar = b2; b3bar = b3;
        % Coupled soft shrinkage update of b1
        s1x = b1x + u(:,[2:end end],:)-u; 
        s1y = b1y + u([2:end end],:,:)-u;
        s1z = b1z + resz*(u(:,:,[2:end end])-u);
        temp = sqrt(s1x.^2 + s1y.^2 + s1z.^2);
        t = sign(temp) .* max(abs(temp)-mu(1)/sigma,0);
        s1x = s1x .* (t ./ max(temp,1e-9));
        s1y = s1y .* (t ./ max(temp,1e-9));
        s1z = s1z .* (t ./ max(temp,1e-9));
        b1x = b1x + u(:,[2:end end],:)-u - s1x;
        b1y = b1y + u([2:end end],:,:)-u - s1y; 
        b1z = b1z + resz.*(u(:,:,[2:end end])-u) - s1z;
        % Soft shrinkage update of b2
        s2 = b2 + s([2:end end],:,:) - s;
        b2 = s2 - sign(s2).* max(abs(s2)-lambda/sigma,0);
        % Soft shrinkage für y3 und b3 - use s2 for storage optimization
        s2 = b3 + s;
        b3 = s2 - sign(s2).* max(abs(s2)-mu(2)/sigma,0);

        % Update Step:
        b1xbar = 2*b1x - b1xbar;
        b1ybar = 2*b1y - b1ybar;
        b1zbar = 2*b1z - b1zbar;
        b2bar = 2*b2 - b2bar;
        b3bar = 2*b3 - b3bar;

        % Increase counter
        k = k+1;
    end
    if verbose
        fprintf('\rStripe removal in %d steps completed.\n',k-1)
    end
    if GPU
        u = gather(u);
    end
end