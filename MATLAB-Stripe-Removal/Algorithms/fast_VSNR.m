%% Variational Stationary Noise Remover (VSNR)
% [Fehrenbach et al., 2012, Variational algorithms to remove stationary
% noise. Application to microscopy imaging.]
% The code is publically available on the homepage of Pierre Weiss
% and was modified to reduce outputs and allow GPU usage.

% function u=VSNR(u0,eps,p,Gabor,alpha,maxit,prec,C)
%
% This function helps removing "structured"
% additive noise. By "structured", we mean that the noise consists of
% convolving a white noise with a given filter.
%
% PD stands for Primal-Dual, as the core of the program is a
% first-order primal-dual algorithm (described in
% "A first-order Primal-Dual Algorithm for convex problems with application to imaging",
% by A. Chambolle and T. Pock).
%
% This function solves (in the sense that the duality gap is less than prec):
% Min_{Lambda} ||nabla u||_{1,eps} + sum_{i=1}^m alpha(i)||Lambda_i||_p_i
% over the constraint that ||lambda||_infty<=C
%
% where :
% - u=u0+sum_{i=1}^m conv(Lambda_i,Gabor(i))
% - ||q||_{1,eps}= sum_{i=1}^n f(|q|) where f(x)=|x| if |x|>eps and
% f(x)=x^2/2eps+eps/2 otherwise.
% - ||.||_p is the standard l^p-norm
%
% IN :
% - u0 : original image (size nx*ny).
% - eps : regularization parameter for TV norm (can be 0).
% - p : either a value in {1,2,Inf} or a vector of size m with values in the previous set.
% - Gabor : a nx*ny*m array containing the shape of m Gabor filters.
% - alpha : either a value in R_+^* either a vector of size m in the with
% all value in R_+^*.
% - maxit : maximum iteration number.
% - prec : a value in R_+^* that specifies the desired precision (typical=1e-2).
% - C    : l-infinite constraint on lambda.
% - GPU  : 1 if GPU is available, 0 for CPU
%
% OUT :
% - u : u0+sum_{i=1}^m conv(Lambda_i,Gabor(i))
% - EstP and EstD : solutions of primal and dual problem.
%
% Any comments: please contact Pierre Weiss, pierre.armand.weiss@gmail.com

function u=fast_VSNR(u0,eps,p,Filter,alpha,maxit,prec,C,GPU)
    % Retrieves informations about the problem type. m is the filters number.
    if ndims(Filter)==2
        m=1;
        [nx,ny]=size(Filter);
    else
        [nx,ny,m]=size(Filter);
    end
    
    % Adjust parameter size to m - if not fitting
    if (m>1)
        if (length(p)~=m)
            p=p(1)*ones(m,1);
        end
        if (length(alpha)~=m)
            alpha=alpha(1)*ones(m,1);
        end
    end
    
    % Primal Variable - initialize
    lambda=zeros(nx,ny,m); 
    lambdabar=lambda;
    % Dual Variables - initialize
    qx=zeros(nx,ny); 
    qy=zeros(nx,ny); 
    
    % Computes Fourier transforms of Filters - Efficiency
    FFilter=zeros(size(Filter));
    for i=1:m
        FFilter(:,:,i)=fft2(Filter(:,:,i));
    end
    
    [gap0,~,~] = DualityGap(u0,lambda,qy,qx,FFilter,eps,p,alpha,C,1);
    
    %% Computation of the largest singular value of A
    d1=zeros(size(u0));
    d1(end,1)=1;d1(1,1)=-1;
    d2=zeros(size(u0));
    d2(1,end)=1;d2(1,1)=-1;
    d1h=fft2(d1);
    d2h=fft2(d2);
    
    H=zeros(size(u0));
    for i=1:m
        H=H+abs(FFilter(:,:,i)).^2;
    end
    L=sqrt(max(H(:).*(abs(d1h(:)).^2+abs(d2h(:)).^2)));
    %disp(sprintf('Operator norm : %1.5f',L))
    clear d1 d2 d1h d2h H;
    
    u=u0;
    gap=Inf;
    
    gamma=min(alpha(:));
    weight=1;
    tau=weight/L;
    sigma=1/(tau*L^2);
    theta=1;
    n=2;
    
    if (isnan(gap0))
        fprintf('INITIAL DUAL GAP IS INFINITE -- PROGRAM WILL STOP AFTER %i ITERATIONS \n',maxit);
        gap0=1e16;
    end
    
    if GPU
        % Transform everything to gpuarray
        lambda = gpuArray(lambda);
        u0 = gpuArray(u0);
        qx = gpuArray(qx);
        qy = gpuArray(qy);
    end

    %% The actual algorithm
    while (n<maxit+1)
        lambdaold = lambda;
        % Step 1: Update q
        b=gpuArray(zeros(size(u0)));
        for i=1:m
            b=b+ifft2(fft2(lambdabar(:,:,i)).*FFilter(:,:,i));
        end
        b=real(b); % at this point b represents the noise.
        
        % Gradient (corresponds to tilde q_n in the article)
        qytilde = qy+sigma*(u0([2:end end],:) - u0) + sigma*(b([2:end end],:) - b);
        qxtilde = qx+sigma*(u0(:,[2:end end]) - u0) + sigma*(b(:,[2:end end]) - b);
        
        % Resolvent operator... - Das stimmt doch nicht mit dem Paper
        % überein
        nqq=sqrt(qytilde.^2+qxtilde.^2);
        qy=qytilde./(max(nqq,eps*sigma+1));
        qx=qxtilde./(max(nqq,eps*sigma+1));
        
        % Step 2: Update lambda
        % D^T q
        dyq = qy([1 1:end-1],:) - qy;
        dyq(1,:) = - qy(1,:);
        dyq(end,:) = qy(end-1,:);
        dxq = qx(:,[1 1:end-1]) - qx;
        dxq(:,1) = - qx(:,1);
        dxq(:,end) = qx(:,end-1);
        % 
        if GPU
            Astarq=gpuArray(zeros(size(lambda)));
        else 
            Astarq = zeros(size(lambda));
        end
        FgradTq = fft2(dxq+dyq);
        for i=1:m
            Astarq(:,:,i)=real(ifft2(conj(FFilter(:,:,i)).*FgradTq));
        end
        lambda=lambda-tau*Astarq;
    
        %Computation of the resolvent of (I+tau partial G)^{-1}
        for i=1:m
            lambda(:,:,i)=Prox_Phi(lambda(:,:,i),p(i),tau,alpha(i),C);
        end
    
      % Update Step size / parameters
        if (sum(p==2)==length(p)) %If all phi_i are l2
            if (eps>0)
                mu=2*sqrt(gamma*eps)/L;
                tau=mu/(2*gamma);
                sigma=mu/(2*eps);
                theta=1/(1+mu);
            else
                theta=1/sqrt(1+2*gamma*tau);
                tau=theta*tau;
                sigma=sigma/theta;
            end
        % else theta=1, tau=tau, sigma = sigma;
        end
    
        % Extrapolation
        lambdabar=lambda+theta*(lambda-lambdaold);
    
        % Display state of optimization
        if (mod(n,10)==0)
            fprintf('\rIteration: %d / %d', n, maxit)
        end
        % Increase counter
        n=n+1;
    end
    
    b=zeros(size(u0));
    for i=1:m
        b=b+ifft2(fft2(lambda(:,:,i)).*FFilter(:,:,i));
    end
    b=real(b); %at this point b represents the noise.
    %Current estimate of the denoised image
    if GPU
        u=gather(u0+b);
    else
        u = u0 + b;
    end

end
%% Companion functions

%This function computes Phi^*(lambda)
%where:
%Phi(x)=||x||_1 if p=1
%Phi(x)=1/2||x||_2^2 if p=2
%Phi(x)=0 if p=Infty and ||x||_inf<=1 inf otherwise
function v=PhiStar(lambda,p,C,alpha)
    if (p==1)
      v=sum(max(0,C*(abs(lambda(:))-alpha)));
    elseif (p==2)
      v=alpha*(sum(min(abs(lambda(:)/alpha),C).*abs(lambda(:)/alpha)- 0.5*min(abs(lambda(:)/alpha),C).^2));
    elseif (p==Inf)
      v=min(C,alpha)*sum(abs(lambda(:)));
    end
end

function n=Phi(x,p,alpha)
    if (p==1)
      n=alpha*sum(abs(x(:)));
    elseif (p==2)
      n=alpha/2*sum(x(:).^2);
    elseif (p==Inf)
      if max(abs(x(:)))>alpha
        n=Inf;
      else
        n=0;
      end
    end
end

%This function solves :
% argmin_{|y|<=C} tau alpha ||y||_p + 1/2 ||M(y-x)||_2^2
function y=Prox_Phi(x,p,tau,alpha,C)
    if (tau==0)
      y=x./(max(1,abs(x)/C));
      return
    end
    if (p==1)
      tau=alpha*tau;
      y=max(abs(x)-tau,0);
      y=y.*sign(x)./(max(abs(y)/C,1));
    elseif (p==2)
      tau=tau*alpha;
      y=x./(tau+1);
      y=y./(max(1,abs(y)/C));
    elseif (p==Inf)
      delta=min(alpha,C);
      y=x./(max(1,abs(x)/delta));
    end
end

% This function calculates the duality gap
function [gap,P,D] = DualityGap(u0,lambda,qy,qx,FFilter,eps,p,alpha,C,n)
    s = zeros(size(u0));
    for k=1:size(FFilter,3)
        s = s+ifft2(fft2(lambda(:,:,k)).*FFilter(:,:,k));
    end
    s=real(s); %at this point s represents the noise.
    
    %Computation of the primal cost in gap
    dy = u0([2:end end],:) - u0 + s([2:end end],:) - s;
    dx = u0(:,[2:end end]) - u0 + s(:,[2:end end]) - s;

    if (eps==0)
        P = sum(sqrt(dy(:).^2 + dx(:).^2));
    else
        P = dy(:).^2 + dx(:).^2;
        P = sum(min(P/eps,sqrt(P)) - .5*min (P/eps,eps));
    end
    for k=1:size(FFilter,3)
        P = P + Phi(lambda(:,:,k),p(k),alpha(k));
    end

    %Computation of the dual cost
    dyT = qy([1 1:end-1],:) - qy;
    dyT(1,:) = -qy(1,:);
    dyT(end,:) = qy(end-1,:);
    dxT = qx(:,[1 1:end-1]) - qx;
    dxT(:,1) = - qx(:,1);
    dxT(:,end) = qx(:,end-1);

    ATq = zeros(size(lambda));
    tmp = fft2(dyT+dxT);
    for k=1:size(FFilter,3)
      ATq(:,:,k)=real(ifft2(conj(FFilter(:,:,k)).*tmp));
    end
    D = 0;
    for k=1:size(FFilter,3)
      D = D - PhiStar(-ATq(:,:,k),p(k),C,alpha(k));
    end
    %Computation of Fstar=F*(q)
    tmp = sqrt(qy.^2+qx.^2);
    if (max(tmp(:))>1)
        if n == 1
            D = -Inf;
        end
        % else D = D - 0;
    else
      D = D - eps/2*norm(tmp(:))^2 ...
            + sum((u0([2:end end],:)-u0).*qy,'all') ...
            + sum((u0(:,[2:end end])-u0).*qx,'all');
    end
    gap = P - D;
end