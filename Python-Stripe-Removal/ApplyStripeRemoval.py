from PIL import Image
import tifffile as tiff
import numpy as np
import torch
import matplotlib.pyplot as plt
from GeneralStripeRemover import GeneralStripeRemover


if __name__ == '__main__':
    # Open 2D image
    # img = Image.open("CircularTestingPattern_Reference.png")
    # img = torch.from_numpy(np.array(img)).float()
    # Open 3D image
    img = tiff.imread("SphericalTestingPattern_256_8.tiff")
    img = torch.from_numpy(img).float()

    # Normalize image
    img = (img - img.min())/(img.max() - img.min())

    # Process image
    result = GeneralStripeRemover(img,iterations=5000,proj=True,mu=[10,0.1],resz=1,direction=[1.,0.,0.],verbose=True)

    # Visualize result
    if img.dim() == 2:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
        ax1.imshow(img,cmap='gray')
        ax2.imshow(result,cmap='gray')
        plt.tight_layout()
        plt.show()
    elif img.dim() == 3:
        fig, axes = plt.subplots(3, 2, figsize=(10, 10))
        axes[0,0].imshow(img[:,:,img.shape[2]//2],cmap='gray')
        axes[0,1].imshow(result[:,:,img.shape[2]//2],cmap='gray')
        axes[1,0].imshow(img[:,img.shape[1]//2,:],cmap='gray')
        axes[1,1].imshow(result[:,img.shape[1]//2,:],cmap='gray')
        axes[2,0].imshow(img[img.shape[0]//2,:,:],cmap='gray')
        axes[2,1].imshow(result[img.shape[0]//2,:,:],cmap='gray')
        plt.tight_layout()
        plt.show()