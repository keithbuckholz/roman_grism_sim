import numpy as np
from astropy.io import fits
import sys
import os

valid_args = ["like_fits", "save_to", "make_seg", "make_direct"]

def make_seg(**kwargs):
    """
    Makes an empty segmentation file. Opens template fits. Copies
    header down. Creates a similar zero array, and save fits.
    """

    # Open fits, store header, create zero data array
    open_fits = fits.open(kwargs["like_fits"])
    hdr = open_fits[1].header
    data = np.zeros_like(open_fits[1].data)

    # Save fits file
    filename = os.path.join(kwargs["save_to"], "empty_seg.fits")
    fits.writeto(filename, data=data, header=hdr, overwrite=True)

    return 0

def make_direct(**kwargs):
    """
    Makes an empty direct image. Opens template fits. Creates a similar zero_array.
    Craft HDUList and save fits. 
    """

    # Open fits, store header, create zero data array
    open_fits = fits.open(kwargs["like_fits"])
    hdr = open_fits[1].header
    data = np.zeros_like(open_fits[1].data)

    # Insert GRIZLI header info
    open_fits[0].header["INSTRUME"] = "ROMAN"
    open_fits[0].header["FILTER"] = "det1"
    hdr["CONFFILE"] = "/Users/keith/astr/research_astr/roman_configs/Roman.det1.07242020.conf"

    # Place zero data and modified SCI HDU header into ImageHDU
    imagehdu = fits.ImageHDU(data=data, header=hdr, name="SCI")

    # Create HDUList
    hdul = [open_fits[0], imagehdu, open_fits[2], open_fits[3]]

    # Save fits file
    filename = os.path.join(kwargs["save_to"], "empty_direct.fits")
    fits.HDUList(hdul).writeto(filename, overwrite=True)

    return 0

def main():
    """
    Process kwargs and call make_seg and/or make_direct
    """

    # Catch sys.argv and set default kwargs
    argv = sys.argv[1:]
    kwargs = {"make_seg": False,
              "make_direct": False}

    # Parse sys.argv
    for arg in argv:
        try:
            # Split key,value pairs
            key, value = arg.split("=")
            if (key == "make_seg") or (key == "make_direct"):
                value = bool(value)

        # Raise Exception if command line argument format is not valid
        except ValueError as e:
            msg = f"{arg} is not valid. Args must be key=value format."
            raise e from ValueError(msg)
        
        # Check that provided command line arguments are valid
        if key not in valid_args:
            raise TypeError(f"Encountered unexpected argument: {key}")
        
        kwargs[key] = value
    
    # Make requested files using their functions
    if kwargs["make_direct"]:
        make_direct(**kwargs)
    
    if kwargs["make_seg"]:
        make_seg(**kwargs)

    return 0

if __name__ == "__main__":
    main()