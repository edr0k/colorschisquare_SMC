from astropy.io import fits


def fz2fits(image):
    """
    It converts SPLUS images
    from .fz to .fits
    """
    datos = fits.open(image)[1].data
    heada = fits.open(image)[1].header
    imageout = image[:-2] + 'fits'
    print
    'Creating file: '
    print
    imageout
    fits.writeto(imageout, datos, heada, clobber=True)


fz2fits('STRIPE82-0119/STRIPE82-0119_F378_swp.fz')
fz2fits('STRIPE82-0119/STRIPE82-0119_F395_swp.fz')
fz2fits('STRIPE82-0119/STRIPE82-0119_F410_swp.fz')
fz2fits('STRIPE82-0119/STRIPE82-0119_F430_swp.fz')
fz2fits('STRIPE82-0119/STRIPE82-0119_F515_swp.fz')
fz2fits('STRIPE82-0119/STRIPE82-0119_F660_swp.fz')
fz2fits('STRIPE82-0119/STRIPE82-0119_F861_swp.fz')
fz2fits('STRIPE82-0119/STRIPE82-0119_G_swp.fz')
fz2fits('STRIPE82-0119/STRIPE82-0119_I_swp.fz')
fz2fits('STRIPE82-0119/STRIPE82-0119_R_swp.fz')
fz2fits('STRIPE82-0119/STRIPE82-0119_U_swp.fz')
fz2fits('STRIPE82-0119/STRIPE82-0119_Z_swp.fz')
