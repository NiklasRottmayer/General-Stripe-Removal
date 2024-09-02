# Paper
https://www.biorxiv.org/content/10.1101/2024.02.02.578531v2


# Summary
This repository contains Matlab code for removing stripe artifacts from 2D or 3D images and python code for generating synthetic light-sheet microscopy images. 

Corrupted Image | Destriping Result
:-------------------------:|:-------------------------:
![grafik](https://github.com/user-attachments/assets/19550398-cffe-455d-928c-1d4d94946fd9)| ![grafik](https://github.com/user-attachments/assets/087da623-b4a6-4c6f-a29f-c043255e7fa2)
![15SnBronze_unsheared_crop_Slice23](https://github.com/user-attachments/assets/2291079b-2e22-4866-ab3f-6df695ae5ed2) | ![GSR3D_mu3-0p05_steps25000_proj1_resz0p33_direction1-0-0_normalize0_15SnBronze_Slice25_0p25-1p1](https://github.com/user-attachments/assets/ea3cdf22-d2f8-484b-82ae-ce7e6449afb0)
![TestImage_t3_l256](https://github.com/user-attachments/assets/6a55ad7c-9402-42be-bb6d-9c9b2329810c) | ![GSR2D_mu7-0p1_steps25000_proj1_resz0_direction1-0_normalize0_TestImage-t3-l256_prep](https://github.com/user-attachments/assets/7c77239d-da7f-4293-9a03-de9d9d3515be)


## Stripe Removal
For removing stripe artifacts we offer the code to four different methods including our own solution which generalizes well to all tested settings. It extends on work by Liu et al. [1]. For the selection of parameters, see the provided scripts and algorithms. 

1. (ours) **General Stripe Remover**: A variational method that penalizes undesired image and stripe features to remove artifacts. The objective function is
   
  $$\underset{u+s=u_0}{\mathrm{argmin}}\quad \mu_1 \lVert\nabla_{x,y,z} u\rVert_{2,1}+ \lVert\nabla_y s\rVert_1 + \mu_2 \lVert s \rVert_1 + \iota_{[0,1]^N}(u) $$
  
  which is optimized using the primal-dual gradient hybrid method with extrapolation of the dual variable (PDHGMp) [2]. It supports the choice of stripe direction in 3D and a variable resolution in stack (z-)direction.  

2. (MDSR*) **Multi-directional Stripe Remover**: A powerful Fourier filtering approach based on the non-subsampled contourlet transform, a multi-scale and directional decomposition. Compared to the original proposition [3] we restrict filtering to subimages of a limited range of directions. Only slice-wise processing is currently supported.

3. (VSNR) **Variational Stationary Noise Remover**: A versatile variational method capable of removing repeating artifacts through the choice of patterns. Only slice-wise processing is currently supported [4].

4. (WFF) **Wavelet-Fourier Filtering**: one of the earliest successfull implementations of Fourier filtering for removing stripes based on the wavelet decomposition. It is commonly referenced and compared with in literature. [5]

[1] X. Liu, X. Lu, H. Shen, et al., “Oblique stripe removal in remote sensing images via oriented variation,” CoRR abs/1809.02043 (2018).

[2] M. Burger, A. Sawatzky, and G. Steidl, “First order algorithms in variational image processing,” in Splitting Methods in Communication, Imaging, Science, and Engineering, R. Glowinski, S. J. Osher, and W. Yin, eds. (Springer International Publishing, Cham, 2016), pp. 345–407.

[3] X. Liang, Y. Zang, D. Dong, et al., “Stripe artifact elimination based on nonsubsampled contourlet transform for light sheet fluorescence microscopy,” J. Biomed. Opt. 21, 1–9 (2016).

[4] J. Fehrenbach, P. Weiss, and C. Lorenzo, “Variational algorithms to remove stationary noise: Applications to microscopy imaging,” IEEE Trans. on Image Process. 21, 4420–4430 (2012).

[5] B. Münch, P. Trtik, F. Marone, and M. Stampanoni, “Stripe and ring artifact removal with combined wavelet — Fourier filtering,” Opt. Express 17, 8567–8591 (2009).

## Synthetic Data Generation
For the generation of synthetic light-sheet (fluorescence) microscopy data we rely on the Biobeam package [6] which allows us to simulate the propagation of light through 3D space given a refractive index distrubtion of the geometry and other parameters such as wavelength of the light, numerical apartures and voxel size. As geometry we combine several techniques including the fixed placement of solids and the random placement of non-overlapping spheres. For the choice of refractive indices we refer to [7]. 

Corrupted Image | Ground Truth
:-------------------------:|:-------------------------:
![Synthetic-Embryo-flat_prep96-96](https://github.com/user-attachments/assets/8e4207ae-e3c3-46c7-9265-d462b496943d) | ![Synthetic-Embryo-flat_Ideal_prep96-96](https://github.com/user-attachments/assets/9c2bbcc4-48aa-44fd-8813-66fa2c4eb20c)



[6] M. Weigert, K. Subramanian, S. T. Bundschuh, et al., “Biobeam—multiplexed wave-optical
simulations of light-sheet microscopy,” PLoS Comput. Biol. 14, 1–11 (2018).

[7] P. Y. Liu, L. K. Chin, W. Ser, et al., “Cell refractive index for cell biology and disease diagnosis:
past, present and future,” Lab Chip 16, 634–644 (2016).
