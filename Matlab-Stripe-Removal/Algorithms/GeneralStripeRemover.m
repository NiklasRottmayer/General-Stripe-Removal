%% General Stripe Remover
%
% The "GeneralStripeRemover" performs the optimization for the variational
% method
%
% argmin mu1*||D_(x,y,z) u||_(2,1) + ||D_y s||_1 + mu2*||s||_1 + iota_[0,1](u).
%  u+s=F
%
% using the primal-dual gradient hybrid method with extrapolation of the
% dual variable (PDHGMp) [Burger et al., 2014, First Order Algorithms in
% Variational Image Processing] It supports the choice of stripe direction 
% in 3D and a variable resolution in stack (z-)direction.  
%
%-------------------------------------------------------------------------
% Author: Niklas Rottmayer
% Date: 29.08.2024
%-------------------------------------------------------------------------
%
% Input: 
% F             -   corrupted image
% iterations    -   number of steps to optimize
% mu            -   weighting parameters (array[2])
% proj          -   project onto [0,1] (true/false)
% resz          -   ratio of resolutions z-axis to x and y (in [0,1]) 
%                   D_z = resz*(F(:,:,i+1) - F(:,:,i))
% normalize     -   normalize the input image to [0,1] (true/false) 
% direction     -   stripe direction (array[2]/array[3]) 
% verbose       -   print process information (true/false)
%
% Output:
% u             -   Result of stripe removal
%-------------------------------------------------------------------------
% Comments: A stopping criterion could be added in the future.

function u = GeneralStripeRemover(F,iterations,mu,proj,resz,normalize,...
                                  direction,verbose)

    %% Default Parameters

    if nargin < 2, iterations = 5000;   end
    if nargin < 3, mu = [10,0.1];       end
    if nargin < 4, proj = 1;            end
    if nargin < 5, resz = 0;            end
    if nargin < 6, direction = [1,0,0]; end
    if nargin < 7, verbose = 1;         end
    GPU = canUseGPU();
    iterations = ceil(iterations);
    
    %% General Sanity checks

    assert(length(size(F)) >= 2,'The input image has invalid dimensions.') % does not catch arrays / numbers!
    assert(iterations >= 0,'The number of iterations must be non-negative.') 
    assert(length(mu) >= 2,'mu must be an array of length 2.')
    assert(all(mu(1:2) > 0),'mu must contain positive values.')
    assert(proj == 0 || proj == 1,'proj must be boolean (true/false or 1/0).')
    assert(resz >= 0 && resz <= 1, 'resz must be in [0,1].')
    assert(~all(direction == 0),'direction cannot be all zeros.')

    %% Specific Sanity Checks

    if ismatrix(F)
        assert(resz == 0,'resz must be 0 for 2D images.')
        assert(direction(1)~= 0 || direction(2)~= 0,'Processing in 2D requires a valid direction')
    end
    if resz > 0 && resz < 1
        assert(direction(2) == 0,'direction in z-direction with resz in (0,1) is not supported.')
    end

    %% Preparation

    if verbose, fprintf('Intializing Stripe Removal\n... please wait ...'), end

    % Preparing 2D Processing
    if resz == 0
        tau = 0.35; 
        direction = direction(1:2)./norm(direction(1:2));

        % Transformations
        if direction(1) < 0, F = flip(F,1); end
        if direction(2) < 0, F = flip(F,2); end

        abs_direction = abs(direction);
        if abs(direction(2)) > abs(direction(1))
            F = permute(F,[2,1,3]);
            abs_direction = flip(abs_direction,2);
        end
        
        % Determine closest direction
        supported_directions = [[1,0];[2,1];[1,1]];
        supported_directions = supported_directions./vecnorm(supported_directions,2,2);
        [~,dir_case] = min(vecnorm(supported_directions - abs_direction,2,2));

    % Preparing 3D Processing
    else
        tau = 0.28;
        direction = direction(1:3)./norm(direction(1:3));

        % Transformations
        if direction(1) < 0, F = flip(F,1); end
        if direction(2) < 0, F = flip(F,2); end
        if direction(3) < 0, F = flip(F,3); end

        abs_direction = abs(direction);
        [abs_direction,I] = sort(abs_direction,'descend');
        F = permute(F,I);  

        % Determine closest direction
        supported_directions = [[1,0,0];[2,1,0];[1,1,0];[1,1,1];[2,1,1];[2,2,1]];
        supported_directions = supported_directions./vecnorm(supported_directions,2,2);
        [~,dir_case] = min(vecnorm(supported_directions - abs_direction,2,2));
    end


    %% Initialization

    [nx,ny,nz] = size(F);
    sigma = tau; lambda = 30; % Comment: lambda rescales optimization functional.

    if normalize, F = rescale(F,0,1); end 
    if GPU, b1x = gpuArray(zeros([nx,ny,nz])); F = gpuArray(F);
    else,   b1x = zeros([nx,ny,nz]); end

    % helper variables
    b1xbar = b1x; b1y = b1x; b1ybar = b1x;
    b2 = b1x; b2bar = b1x; b3 = b1x; b3bar = b1x;
    if resz > 0, b1z = b1x; b1zbar = b1x; end
    k = 1;
    u = F; s = b1x;

    %% Iteration

    while k <= iterations % Space for stopping criterion 
        if verbose, fprintf('\rIteration: %d / %d', k, iterations), end
           
        % Part 1: Update u and s
        % Image: s1x = D_x^T b1xbar  
        s1x = b1xbar(:,[1 1:end-1],:) - b1xbar;
        s1x(:,1,:) = -b1xbar(:,1,:);
        s1x(:,end,:) = b1xbar(:,end-1,:);

        % Image: s1y = D_y^T b1ybar
        s1y = b1ybar([1 1:end-1],:,:) - b1ybar;
        s1y(1,:,:) = -b1ybar(1,:,:);
        s1y(end,:,:) = b1ybar(end-1,:,:);
        if resz > 0
            % s1z = D_z^T b1zbar
            s1z = resz.*(b1zbar(:,:,[1 1:end-1]) - b1zbar);
            s1z(:,:,1) = -resz.*b1zbar(:,:,1);
            s1z(:,:,end) = resz.*b1zbar(:,:,end-1);
        end

        % Stripes: s2 = D_Theta^T b2bar
        switch dir_case
            case 1 % Adjoint 0° (vertical)
                s2 = b2bar([1 1:end-1],:,:) - b2bar;
                s2(1,:,:) = -b2bar(1,:,:);
                s2(end,:,:) = b2bar(end-1,:,:);
            case 2 % Adjoint 26.6°
                s2 = b2bar([1 2 1:end-2],[1 1:end-1],:) - b2bar;
                s2(1:2,:,:) = -b2bar(1:2,:,:);
                s2(end-1:end,2:end,:) = b2bar(end-3:end-2,1:end-1,:);
                s2(:,1,:) = -b2bar(:,1,:);
                s2(3:end,end,:) = b2bar(1:end-2,end-1,:);
                s2(1:2,end,:) = 0;
                s2(end-1:end,1,:) = 0; 
            case 3 % Adjoint 45°
                s2 = b2bar([1 1:end-1],[1 1:end-1],:) - b2bar;
                s2(1,:,:) = -b2bar(1,:,:); 
                s2(end,2:end,:) = b2bar(end-1,1:end-1,:);
                s2(:,1,:) = -b2bar(:,1,:);
                s2(2:end,end,:) = b2bar(1:end-1,end-1,:);
                s2(1,end,:) = 0;
                s2(end,1,:) = 0;
            case 4 % Space diagonal
                s2 = b2bar([1 1:end-1],[1 1:end-1],[1 1:end-1]) - b2bar;
                s2(1,:,:) = -b2bar(1,:,:); 
                s2(end,2:end,2:end) = b2bar(end-1,1:end-1,1:end-1);
                s2(:,1,:) = -b2bar(:,1,:);
                s2(2:end,end,2:end) = b2bar(1:end-1,end-1,1:end-1);
                s2(:,:,1) = -b2bar(:,:,1);
                s2(2:end,2:end,end) = b2bar(1:end-1,1:end-1,end-1);
                s2(end,1,1) = 0;
                s2(1,end,1) = 0;
                s2(1,1,end) = 0;
            case 5 % Space Off-diagonal 1
                s2 = b2bar([1 2 1:end-2],[1 1:end-1],[1 1:end-1]) - b2bar;
                s2(1:2,:,:) = -b2bar(1:2,:,:); 
                s2(end-1:end,2:end,2:end) = b2bar(end-3:end-2,1:end-1,1:end-1);
                s2(:,1,:) = -b2bar(:,1,:);
                s2(3:end,end,2:end) = b2bar(1:end-2,end-1,1:end-1);
                s2(:,:,1) = -b2bar(:,:,1);
                s2(3:end,2:end,end) = b2bar(1:end-2,1:end-1,end-1);
                s2(end-1:end,1,1) = 0;
                s2(1,end,1) = 0;
                s2(1,1,end) = 0;
            case 6 % Space Off-diagonal 2
                s2 = b2bar([1 2 1:end-2],[1 2 1:end-2],[1 1:end-1]) - b2bar;
                s2(1:2,:,:) = -b2bar(1:2,:,:); 
                s2(end-1:end,3:end,2:end) = b2bar(end-3:end-2,1:end-2,1:end-1);
                s2(:,1:2,:) = -b2bar(:,1:2,:);
                s2(3:end,end-1:end,2:end) = b2bar(1:end-2,end-3:end-2,1:end-1);
                s2(:,:,1) = -b2bar(:,:,1);
                s2(3:end,3:end,end) = b2bar(1:end-2,1:end-2,end-1);
                s2(end-1:end,1,1) = 0;
                s2(1,end-1:end,1) = 0;
                s2(1,1,end) = 0;
        end

        if resz == 0, u = u - tau * sigma * (s1x + s1y);
        else, u = u - tau * sigma * (s1x + s1y + s1z); end
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
        b1xbar = b1x; b1ybar = b1y; b2bar = b2; b3bar = b3;
        if resz > 0, b1zbar = b1z; end

        % Coupled soft shrinkage update of b1
        s1x = b1x + u(:,[2:end end],:)-u; 
        s1y = b1y + u([2:end end],:,:)-u;
        if resz > 0, s1z = b1z + resz*(u(:,:,[2:end end])-u); 
            temp = sqrt(s1x.^2 + s1y.^2 + s1z.^2);
        else
            temp = sqrt(s1x.^2 + s1y.^2);
        end
        t = sign(temp) .* max(abs(temp)-mu(1)/sigma,0);
        s1x = s1x .* (t ./ max(temp,1e-9));
        s1y = s1y .* (t ./ max(temp,1e-9));
        if resz > 0, s1z = s1z .* (t ./ max(temp,1e-9)); end
        b1x = b1x + u(:,[2:end end],:)-u - s1x;
        b1y = b1y + u([2:end end],:,:)-u - s1y; 
        if resz > 0, b1z = b1z + resz.*(u(:,:,[2:end end])-u) - s1z; end

        % Soft shrinkage update of b2
        switch dir_case
            case 1 % vertical stripes
                s1x = s([2:end end],:,:) - s;
            case 2 % 26.6°
                s1x = s([3:end end-1 end],[2:end end],:) - s;
                s1x(end-1:end,:,:) = 0;
                s1x(:,end,:) = 0;
            case 3 % 45° 
                s1x = s([2:end end],[2:end end],:) - s;
                s1x(end,:,:) = 0;
                s1x(:,end,:) = 0; 
            case 4
                s1x = s([2:end end],[2:end end],[2:end end]) - s;
                s1x(end,:,:) = 0;
                s1x(:,end,:) = 0;
                s1x(:,:,end) = 0;
            case 5
                s1x = s([3:end end-1 end],[2:end end],[2:end end]) - s;
                s1x(end-1:end,:,:) = 0;
                s1x(:,end,:) = 0;
                s1x(:,:,end) = 0;
            case 6
                s1x = s([3:end end-1 end],[3:end end-1 end],[2:end end]) - s;
                s1x(end-1:end,:,:) = 0;
                s1x(:,end-1:end,:) = 0;
                s1x(:,:,end) = 0;
        end
        s2 = b2 + s1x;
        b2 = s2 - sign(s2).* max(abs(s2)-lambda/sigma,0);

        % Soft shrinkage of b3
        s2 = b3 + s;
        b3 = s2 - sign(s2).* max(abs(s2)-mu(2)/sigma,0);

        % Update Step:
        b1xbar = 2*b1x - b1xbar;
        b1ybar = 2*b1y - b1ybar;
        if resz > 0, b1zbar = 2*b1z - b1zbar; end
        b2bar = 2*b2 - b2bar;
        b3bar = 2*b3 - b3bar;

        % Increase counter
        k = k+1;
    end

    %% Backtransformation

    if resz == 0
        if abs(direction(2)) > abs(direction(1)), u = permute(u,[2,1,3]); end
        if direction(2) < 0, u = flip(u,2); end
        if direction(1) < 0, u = flip(u,1); end
    else    
        Iinv = [find(I==1),find(I==2),find(I==3)];
        u = permute(u,Iinv); 
    
        if direction(3) < 0, u = flip(u,3); end
        if direction(2) < 0, u = flip(u,2); end
        if direction(1) < 0, u = flip(u,1); end
    end
    if verbose, fprintf('\rStripe removal in %d steps completed.\n',k-1), end
    if GPU, u = gather(u); end

end






