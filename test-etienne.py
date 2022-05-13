import numpy as np
from astropy.io import fits
from scipy.ndimage import zoom
from tqdm import tqdm

im = fits.getdata('cds.fits.gz')

pix_bin = 64 # pixel size of the binning box

nbin = 4096//pix_bin

percentile = 50

box = np.zeros([nbin,nbin])
for i in tqdm(range(nbin)):
    for j in range(nbin):
        # -1 sigma of distrubution. Should be good to remove illuminated pixels
        box[i,j] = np.nanpercentile(im[i*pix_bin:i*pix_bin+pix_bin,j*pix_bin:j*pix_bin+pix_bin],percentile)

from matplotlib import pyplot as plt
f,ax = plt.subplots(1,2)
ax[0].imshow(box)
ax[0].set_title('box')
lowf = zoom(box,pix_bin)
ax[1].imshow(lowf)
ax[1].set_title('lowf')
plt.show()
im/=lowf
fits.writeto('test.fits',im, overwrite = True)

mask = np.zeros_like(im)
mask[im<0.5] = 1
mask[im>1.3] = 1
mask[0:4] = 0
mask[0:4] = 0

fits.writeto('mask.fits',mask, overwrite = True)