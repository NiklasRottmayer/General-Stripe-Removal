# Summary
This repository contains Matlab code for removing stripe artifacts from 2D or 3D images and python code for generating synthetic light-sheet microscopy images. 

## Stripe Removal
For removing stripe artifacts we offer the code to four different methods including our own solution which generalizes well to all tested settings. It extends on work by Liu et al. [1]. For the selection of parameters, see the provided scripts and algorithms. 

1. (ours) **General Stripe Remover**: A variational method that penalizes undesired image and stripe features to remove artifacts. The objective function is
   
  $$\underset{u+s=u_0}{\mathrm{argmin}}\, \mu_1 \lVert\nabla_{x,y,z} u\rVert_{2,1}+ \lVert\nabla_y s\rVert_1 + \mu_2 \lVert s \rVert_1 + \iota_{[0,1]^N}(u) $$
  
  which is optimized using the primal-dual gradient hybrid method with extrapolation of the dual variable (PDHGMp) [2]. It supports the choice of stripe direction in 3D and a variable resolution in stack (z-)direction.  

3. (MDSR*) **Multi-directional Stripe Remover**: A powerful Fourier filtering approach based on the non-subsampled contourlet transform, a multi-scale and directional decomposition. Compared to the original proposition [3] we restrict filtering to subimages of a limited range of directions. Only slice-wise processing is currently supported.

4. (VSNR) **Variational Stationary Noise Remover**: A versatile variational method capable of removing repeating artifacts through the choice of patterns. Only slice-wise processing is currently supported [4].



[1] X. Liu, X. Lu, H. Shen, et al., “Oblique stripe removal in remote sensing images via oriented variation,” CoRR abs/1809.02043 (2018).

[2] M. Burger, A. Sawatzky, and G. Steidl, “First order algorithms in variational image processing,” in Splitting Methods in Communication, Imaging, Science, and Engineering, R. Glowinski, S. J. Osher, and W. Yin, eds. (Springer International Publishing, Cham, 2016), pp. 345–407.

[3] X. Liang, Y. Zang, D. Dong, et al., “Stripe artifact elimination based on nonsubsampled contourlet transform for light sheet fluorescence microscopy,” J. Biomed. Opt. 21, 1–9 (2016).

[4] J. Fehrenbach, P. Weiss, and C. Lorenzo, “Variational algorithms to remove stationary noise: Applications to microscopy imaging,” IEEE Trans. on Image Process. 21, 4420–4430 (2012).

## Synthetic Data Generation
For the generation of synthetic light-sheet (fluorescence) microscopy data we rely on the Biobeam package [5] which allows us to simulate the propagation of light through 3D space given a refractive index distrubtion of the geometry and other parameters such as wavelength of the light, numerical apartures and voxel size. As geometry we combine several techniques including the fixed placement of solids and the random placement of non-overlapping spheres. For the choice of refractive indices we refer to [6]. 

[5] M. Weigert, K. Subramanian, S. T. Bundschuh, et al., “Biobeam—multiplexed wave-optical
simulations of light-sheet microscopy,” PLoS Comput. Biol. 14, 1–11 (2018).

[6] P. Y. Liu, L. K. Chin, W. Ser, et al., “Cell refractive index for cell biology and disease diagnosis:
past, present and future,” Lab Chip 16, 634–644 (2016).
