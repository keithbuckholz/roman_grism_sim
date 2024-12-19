# General imports
import numpy as np
from scipy import constants
from astropy.io import fits
from astropy.table import Table, join
from grizli.model import GrismFLT
import pysynphot as S
import webbpsf
import webbpsf.roman
import argparse
import os

def spectrum_lookup(stellar_type=None) -> Table:

    # TODO the thing

    spectrum_1D = None
    return spectrum_1D

def chunk(spec, start_wave, end_wave):
    """
    Returns a chunk of the spectrum array between start and end wavelengths.

    Parameters
    ----------
    spec: astropy.table.Table
        Spectrum the chunk should come from.

    start_wave: int

    end_wave: int
    """

    # TODO Fix this algorithm to actually lookup the start and 
    # TODO end_indicies for generalization purposes
    start_index = int((start_wave - 10000) / 5)
    end_index = int((end_wave - 10000) / 5) + 1

    chunk_spec = np.asarray(spec[start_index:end_index])

    return chunk_spec

def dispersion_model(roman, wfi, src, spectrum_overlap, npsfs) -> GrismFLT:
    """
    Produce dipsered model. Saves model in the GrismFLT object. Returns that object.

    Parameters
    ----------

    roman: grizli.model.GrismFLT
        GrismFLT object. The pieces of the dispersion will be saved in the GrismFLT object

    spectrum_overlap: int
        The number of data points that overlap as spectrum segments roll-on/off.
    
    npsfs: int
        The number of distinct PSFs to be used. Also, the number of spectrum segments.
    """

    spectrum_overlap = int(spectrum_overlap * 2) # overlap extent; data points
    npsfs = int(npsfs) # number of distinct psfs

    window_x = np.linspace(0, np.pi, spectrum_overlap)
    front_y = (1-np.cos(window_x)) / 2
    back_y = 1 - front_y

    bins = np.linspace(10000, 20000, npsfs + 1)

    piecemeal_sim = np.zeros((4288,4288))

    input_catcher = []
    
    for ii, start_wave in enumerate(bins[:-1]):

        psf = wfi.calc_psf(monochromatic=(start_wave * (10**-10)), fov_pixels=182, oversample=1, source=src)[0].data
        half_psf_thumb = int(psf.shape[0] / 2) # Used for indexing; center_pixel plus/minus half_psf_thumb(nail)

        direct = np.zeros((4288, 4288))
        direct[(2144-half_psf_thumb): (2144+half_psf_thumb), (2144-half_psf_thumb):(2144+half_psf_thumb)] = psf

        roman.direct.data["SCI"] = direct.astype("float32")
        roman.seg = np.where(roman.direct.data["SCI"], 1, 0).astype("float32")

        end_wave = bins[ii+1]

        start_wave -= (spectrum_overlap * 5) * 0.5 # times 5 from indicies to wavelength
        end_wave += ((spectrum_overlap - 1) * 5) * 0.5

        if start_wave < 10000:
            start_wave = 10000
        
        if end_wave > 20000:
            end_wave = 20000

        chunk_spec = chunk(start_wave, end_wave) # extract relevant part of spectrum
        wave = chunk_spec["col1"]
        flux = chunk_spec["col2"]

        # apodize
        if start_wave != 10000:
            flux[:spectrum_overlap] *= front_y

        if end_wave != 20000:
            flux[-spectrum_overlap:] *= back_y

        roman.compute_model_orders(id=1, mag=1, compute_size=False, size=77, is_cgs=True, store=False, 
                                   in_place=True, spectrum_1d=[wave, flux])

        del chunk_spec
        del flux
        del wave

    return roman

def main(empty_fits_dir=None, spectrum_file=None, bandpass_file=None, magnitude=6,
         spectrum_overlap=10, npsfs=20) -> None:
    """
    "If name==main" parses key words. Then calls main.

    If not provided, the main function calls the spectrum lookup function 
    which pulls a star spectrum template from the _______ library. Then, it
    prepares and calls the dispersion_model function. This function produces
    a grism disperion model using a wavelength-dependent PSF. Finally, main 
    saves the model as a fits file, which can later be accessed later to avoid
    needlessly recomputing a given stellar types model. i.e. This trades 
    compute time for storage space.
    """

    # Read in bandpass array; Will use this to renorm the spectrum
    if bandpass_file:
        df = Table.read(bandpass_file, format='fits')
        bp = S.ArrayBandpass(df["WAVELENGTH"], df["THROUGHPUT"])

    # Read in the spectrum
    if spectrum_file:
        spec = Table.read(spectrum_file, format="ascii")
        src = S.ArraySpectrum(wave=spec["col1"], flux=spec["col2"], waveunits="angstroms", fluxunits="flam")
        
        # If bandpass was provided, here is where we renorm the spectrum; skip renorm if bp does not exist
        if bandpass_file:
            src = src.renorm(magnitude, "abmag", bp)
            src.convert("flam")

    # If spectrum was not provided, lookup a spectrum matching the stellar_type in library
    else:
        # TODO Implement lookup function
        None

    if empty_fits_dir:
        empty_direct = os.path.join(empty_fits_dir, "empty_direct.fits")
        empty_seg = os.path.join(empty_fits_dir, "empty_seg.fits")

    else:
        # TODO Call make_empty.py
        None

    # Instantiate GrismFLT; Fill in with empty direct image and segmentation map
    # Pad is hardcoded for now; I have never needed a number other than 100
    roman = GrismFLT(direct_file=empty_direct, seg_file=empty_seg, pad=100)

    # Instantiate WebbPSF object
    wfi = webbpsf.roman.WFI()

    # Produce the wavelength-depedent PSF Grism Dispersion and save the model as a FITS file
    roman = dispersion_model(roman, wfi, src, spectrum_overlap, npsfs)
    roman.save_model()

    return None

if __name__ == "__main__":
    # Initialize arg parser to assist with command line arguments
    parser = argparse.ArgumentParser()

    parser.add_argument("--empty_fits_dir", type=str, default=None,
                        help="An optional filepath pointing to a directory with empty direct and segmentation fits files.")
    
    parser.add_argument("--spectrum_file", type=str, default=None,
                        help="An optional filepath to the spectrum file.")

    parser.add_argument("--bandpass_file", type=str, default=None,
                        help="An optional filepath to a bandpass file.")

    parser.add_argument("--stellar_type", type=str, default=None,
                        help="If spectrum_file is not provided, a required argument indicating stellar type to lookup.")
                    
    parser.add_argument("--magnitude", type=float, default=6,
                        help="An optional argument: Magnitude of the star in the bandpass provided.")

    parser.add_argument("--spectrum_overlap", type=int, default=10,
                        help="An optional argument: number of data points in the overlapping regions of the spectrum segments.")

    parser.add_argument("--npsfs", type=int, default=20,
                        help="The number of distinct PSFs to be used. Also, the number of spectrum segments.")

    args = parser.parse_args()

    main(empty_fits_dir=args.empty_fits_dir,
         spectrum_file=args.spectrum_file, 
         bandpass_file=args.bandpass_file, 
         magnitude=args.magnitude,
         spectrum_overlap=args.spectrum_overlap, 
         npsfs=args.npsfs)