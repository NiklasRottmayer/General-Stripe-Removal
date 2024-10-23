# GENERAL STRIPE REMOVER
#
# The "GeneralStripeRemover" performs the optimization for the variational method
#
#   argmin mu1*||D_(x,y,z) u||_(2,1) + ||D_y s||_1 + mu2*||s||_1 + iota_[0,1](u).
#    u+s=F
#
# using the primal-dual gradient hybrid method with extrapolation of the dual variable (PDHGMp)
# [Burger et al., 2014, First Order Algorithms in Variational Image Processing].
#
# ----------------------------------------------------------------------------------------------------------------------
# Author: Niklas Rottmayer
# Date: 22.10.2024
# ----------------------------------------------------------------------------------------------------------------------
#
# Input:
# F             -   corrupted image
# iterations    -   number of steps to optimize
# mu            -   weighting parameters (array[2])
# proj          -   project onto [0,1] (true/false)
# resz          -   ratio of resolutions z-axis to x and y (in [0,1])
#                   D_z = resz*(F(:,:,i+1) - F(:,:,i))
# normalize     -   normalize the input image to [0,1] (true/false)
# direction     -   direction of stripes (np.array[2] or np.array[3])
# GPU           -   use of GPU if available (true/false)
# verbose       -   print process information (true/false)
#
# Output:
# u             -   Destriping result
#
# ----------------------------------------------------------------------------------------------------------------------
# Comments: A stopping criterion could be added in the future.
#           If you encounter a problem or error, contact me via niklas.rottmayer@rptu.de or GitHub.

import numpy as np
import torch

def GeneralStripeRemover(F, iterations=5000, mu=[10, 0.1], proj=True, resz=0, normalize=False,
                            direction=[1.,0.,0.], GPU=True, verbose=True):

    # General sanity checks
    assert 2 <= F.dim() <= 3,                   'The input image has invalid dimensions.'
    assert iterations >= 0,                     'The number of iterations must be non-negative.'
    assert len(mu) >= 2,                        'mu must be an array of length 2.'
    assert all(m > 0 for m in mu[:2]),          'mu must contain positive values.'
    assert isinstance(proj, bool),              'proj must be boolean (true/false).'
    assert 0 <= resz <= 1,                      'resz must be in [0,1].'
    assert isinstance(normalize, bool),         'normalize must be boolean (true/false).'
    assert not all(d == 0 for d in direction),  'direction must be non-zero.'
    assert isinstance(GPU, bool),               'GPU must be boolean (true/false).'
    assert isinstance(verbose, bool),           'verbose must be boolean (true/false).'

    # Specific sanity checks
    if F.dim() == 2:
        assert resz == 0,                       'resz must be zero for 2D images.'
        F = F.unsqueeze(2)  # transform image to 3D tensor
    elif F.dim() == 3 and direction[2] != 0:
        assert resz == 1,          'stripe direction in z-direction is not supported with resz in [0,1).'
    # Divide by zero prevention
    prec = 1e-9

    with torch.no_grad():
        if GPU:
            GPU = torch.cuda.is_available()

        # Preparation
        if verbose:
            print('Initializing Stripe Removal\n... please wait ...')

        # Processing in 2D / Slice-by-slice
        if resz == 0:
            tau = 0.35
            direction = direction[:2] / np.linalg.norm(direction[:2])
            # Transformation
            if direction[0] < 0:
                F = F.flip(dims=[0])
            if direction[1] < 0:
                F = F.flip(dims=[1])
            abs_direction = np.abs(direction)
            if abs_direction[1] > abs_direction[0]:
                F = F.permute(1,0,2)
                abs_direction = np.flip(abs_direction)
            # Determine closest supported direction
            supported_directions = np.array([[1.,0.],[2.,1.],[1.,1.]])
            supported_directions /= np.linalg.norm(supported_directions,axis=1,ord=2)[:,np.newaxis]
            dir_case = np.argmin(np.linalg.norm(supported_directions - abs_direction[np.newaxis,:],axis=1))
        # Processing in 3D / full stack
        else:
            tau = 0.28
            direction = direction[:3] / np.linalg.norm(direction[:3])
            # Transformation
            if direction[0] < 0:
                F = F.flip(dims=[0])
            if direction[1] < 0:
                F = F.flip(dims=[1])
            if direction[2] < 0:
                F = F.flip(dims=[2])
            abs_direction = np.abs(direction)
            index = np.argsort(abs_direction)[::-1]
            abs_direction = np.array(abs_direction)[index]
            F = F.permute(tuple(index))
            # Determine closest supported direction
            supported_directions = np.array([[1.,0.,0.],[2.,1.,0.],[1.,1.,0.],[1.,1.,1.],[2.,1.,1.],[2.,2.,1.]])
            supported_directions /= np.linalg.norm(supported_directions,axis=1,ord=2)[:,np.newaxis]
            dir_case = np.argmin(np.linalg.norm(supported_directions - abs_direction[np.newaxis,:],axis=1))

        # Initialization
        sigma = tau
        lmb = 30.
        nx, ny, nz = F.size()

        if normalize:
            F = (F - F.min()) / (F.max() - F.min())
            print('Normalization applied')

        if GPU:
            print('GPU utilized')
            torch.set_default_device('cuda')
            F = F.to('cuda')

        # Helper variables
        b1x = torch.zeros((nx, ny, nz))
        b1xbar = torch.zeros((nx, ny, nz))
        b1y = torch.zeros((nx, ny, nz))
        b1ybar = torch.zeros((nx, ny, nz))
        if resz > 0:
            b1z = torch.zeros((nx, ny, nz))
            b1zbar = torch.zeros((nx, ny, nz))
        b2 = torch.zeros((nx, ny, nz))
        b2bar = torch.zeros((nx, ny, nz))
        b3 = torch.zeros((nx, ny, nz))
        b3bar = torch.zeros((nx, ny, nz))
        if dir_case != 0:
            s2 = torch.zeros((nx, ny, nz))

        u = F.clone().reshape((nx,ny,nz))
        s = torch.zeros((nx, ny, nz))

        #Test = torch.zeros((iterations,nx*ny))

        for k in range(iterations):
            if verbose:
                print(f'\rIteration: {1+k} / {iterations}', end='')

            # Part 1: Update u and s
            s1x = -b1xbar.diff(dim=0,prepend=b1xbar[0,:,:].reshape((1,ny,nz)))
            s1x[0,:,:] = -b1xbar[0,:,:]
            s1x[-1, :] = b1xbar[-2, :]

            s1y = -b1ybar.diff(dim=1,prepend=b1ybar[:,0,:].reshape((nx,1,nz)))
            s1y[:,0,:] = -b1ybar[:,0,:]
            s1y[:,-1,:] = b1ybar[:,-2,:]

            if resz > 0:
                s1z = resz * -b1zbar.diff(dim=2,prepend=b1zbar[:,:,0].reshape((nx,ny,1)))
                s1z[:,:,0] = -resz * b1zbar[:,:,0]
                s1z[:,:,-1] = resz * b1zbar[:,:,-2]

            # Stripes: s2 = D_Theta^T b2bar (vertical)
            if dir_case == 0:  # Adjoint 0° (vertical)
                s2 = -b2bar.diff(dim=0,prepend=b2bar[0,:,:].reshape((1,ny,nz)))
                s2[0,:,:] = -b2bar[0,:,:]
                s2[-1,:,:] = b2bar[-2,:,:]
            elif dir_case == 1:  # Adjoint 26.6°
                s2[2:,1:,:] = b2bar[:-2,:-1,:] - b2bar[2:,1:,:]
                s2[:2,:,:] = -b2bar[:2,:,:]
                s2[-2:,1:,:] = b2bar[-4:-2,:-1,:]
                s2[:,0,:] = -b2bar[:,0,:]
                s2[2:,-1,:] = b2bar[0:-2,-2,:]
                s2[:2,-1,:] = 0
                s2[-2:,0,:] = 0
            elif dir_case == 2:  # Adjoint 45°
                s2[1:,1:,:] = b2bar[:-1,:-1,:] - b2bar[1:,1:,:]
                s2[0,:,:] = -b2bar[0,:,:]
                s2[-1, 1:,:] = b2bar[-2,:-1,:]
                s2[:,0,:] = -b2bar[:,0,:]
                s2[1:,-1,:] = b2bar[:-1,-2,:]
                s2[0,-1,:] = 0
                s2[-1,0,:] = 0
            elif dir_case == 3:  # Space diagonal
                s2[1:,1:,1:] = b2bar[:-1,:-1,:-1] - b2bar[1:,1:,1:]
                s2[0,:,:] = -b2bar[0,:,:]
                s2[-1,1:,1:] = b2bar[-2,:-1,:-1]
                s2[:,0,:] = -b2bar[:,0,:]
                s2[1:,-1,1:] = b2bar[:-1,-2,:-1]
                s2[:,:,0] = -b2bar[:,:,0]
                s2[1:,1:,-1] = b2bar[:-1,:-1,-2]
                s2[-1,0,0] = 0
                s2[0,-1,0] = 0
                s2[0,0,-1] = 0
            elif dir_case == 4:  # Space Off-diagonal 1
                s2[2:,1:,1:] = b2bar[:-2,:-1,:-1] - b2bar[2:,1:,1:]
                s2[:2,:,:] = -b2bar[:2,:,:]
                s2[-2:,1:,1:] = b2bar[-4:-2,:-1,:-1]
                s2[:,0,:] = -b2bar[:,0,:]
                s2[2:,-1,1:] = b2bar[:-2,-2,:-1]
                s2[:,:,0] = -b2bar[:,:,0]
                s2[2:,1:,-1] = b2bar[:-2,:-1,-2]
                s2[-2:,0,0] = 0
                s2[0,-1,0] = 0
                s2[0,0,-1] = 0
            elif dir_case == 5:  # Space Off-diagonal 2
                s2[2:,2:,1:] = b2bar[:-2,:-2,:-1] - b2bar[2:,2:,1:]
                s2[:2,:,:] = -b2bar[:2,:,:]
                s2[-2:,2:,1:] = b2bar[-4:-2,:-2,:-1]
                s2[:,:2,:] = -b2bar[:,:2,:]
                s2[2:,-2:,1:] = b2bar[:-2,-4:-2,:-1]
                s2[:,:,0] = -b2bar[:,:,0]
                s2[2:,2:,-1] = b2bar[:-2,:-2,-2]
                s2[-2:,0,0] = 0
                s2[0,-2:,0] = 0
                s2[0,0,-1] = 0
            else:
                raise ValueError('No case for direction was found. Please check direction.')

            # Compute u
            if resz > 0:
                u -= tau * sigma * (s1x + s1y + s1z)
            else:
                u -= tau * sigma * (s1x + s1y)
            s -= tau * sigma * (s2 + b3bar)

            # Reprojection onto u+s=F
            tmp = F - s - u

            u += 0.5*tmp
            s += 0.5*tmp
            if proj:
                s += (u < 0) * u + (u > 1) * (u-1)
                u = u.clamp(min=0,max=1)

            # Updating helper variables
            b1xbar = b1x.clone()
            b1ybar = b1y.clone()
            if resz > 0:
                b1zbar = b1z.clone()
            b2bar = b2.clone()
            b3bar = b3.clone()

            # Coupled soft-shrinkage
            s1x = b1x + u.diff(dim=0,append=u[-1,:,:].reshape((1,ny,nz)))
            s1y = b1y + u.diff(dim=1,append=u[:,-1,:].reshape((nx,1,nz)))
            if resz > 0:
                s1z = b1z + u.diff(dim=2,append=u[:,:,-1].reshape((nx,ny,1)))
                tmp = (s1x**2 + s1y**2 + s1z**2).sqrt()
            else:
                tmp = (s1x**2 + s1y**2).sqrt()
            t = tmp.sign() * (tmp.abs()-mu[0]/sigma).clamp(min=0)
            s1x *= t/(tmp.clamp(min=prec))
            s1y *= t/(tmp.clamp(min=prec))
            if resz > 0:
                s1z *= t/(tmp.clamp(min=prec))
                b1z += u.diff(dim=2,append=u[:,:,-1].reshape((nx,ny,1))) - s1z
            b1x += u.diff(dim=0,append=u[-1,:,:].reshape((1,ny,nz))) - s1x
            b1y += u.diff(dim=1,append=u[:,-1,:].reshape((nx,1,nz))) - s1y

            # Soft shrinkage of b2
            if dir_case == 0:
                s1x = s.diff(dim=0,append=s[-1,:,:].reshape((1,ny,nz)))
            elif dir_case == 1:
                s1x[:-2,:-1,:] = s[2:,1:,:] - s[:-2,:-1,:]
                s1x[-2:,:,:] = 0
                s1x[:,-1:,:] = 0
            elif dir_case == 2:
                s1x[:-1,:-1,:] = s[1:,1:,:] - s[:-1,:-1,:]
                s1x[-1:,:,:] = 0
                s1x[:,-1:,:] = 0
            elif dir_case == 3:
                s1x[:-1,:-1,:-1] = s[1:,1:,1:] - s[:-1,:-1,:-1]
                s1x[-1:,:,:] = 0
                s1x[:,-1:,:] = 0
                s1x[:,:,-1:] = 0
            elif dir_case == 4:
                s1x[:-2,:-1,:-1] = s[2:,1:,1:] - s[:-2,:-1,:-1]
                s1x[-2:,:,:] = 0
                s1x[:,-1:,:] = 0
                s1x[:,:,-1:] = 0
            elif dir_case == 5:
                s1x[:-2,:-2,:-1] = s[2:,2:,1:] - s[:-2,:-2,:-1]
                s1x[-2:,:,:] = 0
                s1x[:,-2:,:] = 0
                s1x[:,:,-1:] = 0
            else:
                raise ValueError('No case for direction was found. Please check direction.')

            s2 = b2 + s1x
            b2 = s2 - s2.sign() * (s2.abs()-lmb/sigma).clamp(min=0)

            # Soft shrinkage of b3
            s2 = b3 + s
            b3 = s2 - s2.sign() * (s2.abs()-mu[1]/sigma).clamp(min=0)

            # Update step
            b1xbar = 2*b1x - b1xbar
            b1ybar = 2*b1y - b1ybar
            if resz > 0:
                b1zbar = 2*b1z - b1zbar
            b2bar = 2*b2 - b2bar
            b3bar = 2*b3 - b3bar

        # Backtransformation
        if resz == 0:
            if np.abs(direction)[1] > np.abs(direction)[0]:
                u = u.permute(1,0,2)
            if direction[1] < 0:
                u = u.flip(dims=[1])
            if direction[0] < 0:
                u = u.flip(dims=[0])
        else:
            Iinv = np.concatenate([np.where(index == x)[0] for x in range(3)])
            u = u.permute(tuple(Iinv))
            if direction[2] < 0:
                u = u.flip(dims=[2])
            if direction[1] < 0:
                u = u.flip(dims=[1])
            if direction[0] < 0:
                u = u.flip(dims=[0])

    return u.squeeze().cpu()