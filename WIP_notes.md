# WIP Issues

1) Factor of 2 issue; matching scale/count rate 
    - Check scale issue with template; Sum is higher flux than template, max is close to appropriate flux
2) Wavelength-depedent PSF 

0) Producing direct image? (Minor capability for testing; we don't want to scale this)

# Going forward plans

1) [x] Track down factor of 2 error.

2) [x] Check brightest and faintest stars to correct star radius function.

3) [p] Check scale of extracted spectra.

4) [p] Apodize wavelength cuts to roll on/off when attempting wavelength-dependent PSF.

5) [x] Where does Grizli get total flux from? Direct Image? Provided Spectrum? Try doubling direct image values and comparing grism sim scale. Does it change? Do it again with provided spectrum.


# Resolutions

1) Units error. Fixing the error assuages the issue. Switching to use the sed from HLSS products improves things further. Also, padding inclusions offsets seg map locations; need to be sure to track padding presence intentionally.

2) Using the star_radius function as found using logger pro. Wrapped in min(0, max(func, 40)) to ensure a max radius of 40 and min raiuds of 8. Testing in sim_and_extraction supports this solution.

5) Built single_object_playground. Grizli seems to get total flux from spectrum. Disperses total flux according to direct image.

# Progress

3) See single_pix.ipynb. Extracted spectrum from a single pixel object matches expectation (almost perfectly). Expectation is define: (spectrum_template * sensitivity_curve) -> summed in ~10 angstrom bins; binned array contains flux expectations. Also see sim_and_extraction: simulated each object individually to avoid contamination, and extracting in the sim for-loop. Things seem about right. Object 17127 has equal sum and max extractions, but does not sit exactly on the template. Candidate for scale checks.

4) Apodizing makes a big difference. The window needs to be refined. Try larger bins; does fewer discontinuities lower the amplitude of the erroneous wiggles? Using only 5 bins largely reduces the FFT wiggles issues, in exchange for a steps isue. Need to figure out how to smooth. Interpolation? Convolution?

6) Should build simple extraction methods next.

# Quick Notes Space