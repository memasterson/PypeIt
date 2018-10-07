""" Module for LRIS specific codes
"""
from __future__ import absolute_import, division, print_function


import glob

import numpy as np
from astropy.io import fits

from pypeit import msgs
from pypeit.spectrographs import spectrograph
from ..par.pypeitpar import DetectorPar
from pypeit.par.pypeitpar import CalibrationsPar
from .. import telescopes
from pypeit.core import framematch
from pypeit.core import parse
from pypeit import debugger
from pypeit.par import pypeitpar

class GeminiGMOSSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Gemini/GMOS specific code
    """

    def __init__(self):

        # Get it started
        super(GeminiGMOSSpectrograph, self).__init__()
        self.timeunit = 'isot'  # Synthesizes date+time

    def gemini_header_keys(self):
        def_keys = self.default_header_keys()
        def_keys[0]['time'] = 'OBSEPOCH'      # The time stamp of the observation (i.e. decimal MJD)
        def_keys[0]['dispname'] = 'GRATING'      # The time stamp of the observation (i.e. decimal MJD)
        def_keys[0]['idname'] = 'OBSTYPE'     # Frame type
        def_keys[0]['decker'] = 'MASKNAME'
        def_keys[0]['wavecen'] = 'CENTWAVE'
        def_keys[0]['exptime'] = 'EXPTIME'
        def_keys[0]['date'] = 'DATE-OBS'
        def_keys[0]['time'] = 'TIME-OBS'
        def_keys[0]['target'] = 'OBJECT'
        
        def_keys[1] = {}
        def_keys[1]['binning'] = 'CCDSUM'
        # Return
        return def_keys

    def _set_calib_par(self, user_supplied=None):
        self.calib_par = CalibrationsPar()

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype == 'science':
            return good_exp & (fitstbl['idname'] == 'OBJECT')
        if ftype == 'arc':
            return good_exp & (fitstbl['idname'] == 'ARC')
        if ftype == 'pixelflat' or ftype == 'trace':
            return good_exp & (fitstbl['idname'] == 'FLAT')
        if ftype == 'bias':
            return good_exp & (fitstbl['idname'] == 'BIAS')

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

    @staticmethod
    def default_pypeit_par():
        """
        Set default parameters for Keck LRISb reductions.
        """
        par = pypeitpar.PypeItPar()
        # Set wave tilts order
        par['calibrations']['slits']['sigdetect'] = 20.
        par['calibrations']['slits']['polyorder'] = 3
        par['calibrations']['slits']['fracignore'] = 0.02
        par['calibrations']['slits']['pcapar'] = [3,2,1,0]

        # 1D wavelength solution
        par['calibrations']['wavelengths']['rms_threshold'] = 0.40  # Might be grating dependent..

        # Overscan subtract the images
        #par['calibrations']['biasframe']['useframe'] = 'overscan'

        # Alter the method used to combine pixel flats
        par['calibrations']['pixelflatframe']['process']['combine'] = 'median'
        par['calibrations']['pixelflatframe']['process']['sig_lohi'] = [10.,10.]

        # Scienceimage default parameters
        par['scienceimage'] = pypeitpar.ScienceImagePar()
        # Always flux calibrate, starting with default parameters
        par['fluxcalib'] = pypeitpar.FluxCalibrationPar()
        # Always correct for flexure, starting with default parameters
        par['flexure'] = pypeitpar.FlexurePar()

        # Set the default exposure time ranges for the frame typing
        #par['scienceframe']['exprng'] = [30, None]

        return par

    def load_raw_img_head(self, raw_file, det=None, **null_kwargs):
        """
        Wrapper to the raw image reader for LRIS

        Args:
            raw_file:  str, filename
            det: int, REQUIRED
              Desired detector
            **null_kwargs:
              Captured and never used

        Returns:
            raw_img: ndarray
              Raw image;  likely unsigned int
            head0: Header

        """
        raw_img, head0, _ = read_gmos(raw_file, det=det)

        return raw_img, head0


    def get_image_shape(self, filename=None, det=None, **null_kwargs):
        """
        Overrides :class:`Spectrograph.get_image_shape` for LRIS images.

        Must always provide a file.
        """
        # Cannot be determined without file
        if filename is None:
            raise ValueError('Must provide a file to determine the shape of an LRIS image.')

        # Use a file
        self._check_detector()
        self.naxis = (self.load_raw_frame(filename, det=det)[0]).shape
        return self.naxis

    def get_image_section(self, filename, det, section='datasec'):
        """
        Return a string representation of a slice defining a section of
        the detector image.

        Overwrites base class function to use :func:`read_deimos` to get
        the image sections.

        .. todo::
            - It feels really ineffiecient to just get the image section
              using the full :func:`read_deimos`.  Can we parse that
              function into something that can give you the image
              section directly?

        This is done separately for the data section and the overscan
        section in case one is defined as a header keyword and the other
        is defined directly.

        Args:
            filename (str):
                data filename
            det (int):
                Detector number
            section (:obj:`str`, optional):
                The section to return.  Should be either datasec or
                oscansec, according to the :class:`DetectorPar`
                keywords.

        Returns:
            list, bool: A list of string representations for the image
            sections, one string per amplifier, followed by three
            booleans: if the slices are one indexed, if the slices
            should include the last pixel, and if the slice should have
            their order transposed.
        """
        # Read the file
        temp, head0, secs = read_gmos(filename, det=det)
        if section == 'datasec':
            return secs[0], False, False, False
        elif section == 'oscansec':
            # Need to flip these
            return secs[1], False, False, False
        else:
            raise ValueError('Unrecognized keyword: {0}'.format(section))

    def gemini_get_match_criteria(self):
        match_criteria = {}
        for key in framematch.FrameTypeBitMask().keys():
            match_criteria[key] = {}
        # Science
        match_criteria['science']['number'] = 1
        # Standard
        match_criteria['standard']['number'] = 1  # Can be over-ruled by flux calibrate = False
        match_criteria['standard']['match'] = {}
        match_criteria['standard']['match']['decker'] = ''
        # Bias
        match_criteria['bias']['number'] = 5
        match_criteria['bias']['match'] = {}
        match_criteria['bias']['match']['decker'] = ''
        # Pixelflat
        match_criteria['pixelflat']['number'] = 1
        match_criteria['pixelflat']['match'] = {}
        match_criteria['pixelflat']['match']['decker'] = ''
        # Traceflat
        match_criteria['trace']['number'] = 1
        match_criteria['trace']['match'] = {}
        match_criteria['trace']['match']['decker'] = ''
        # Arc
        match_criteria['arc']['number'] = 1
        match_criteria['arc']['match'] = {}
        match_criteria['arc']['match']['decker'] = ''

        # Return
        return match_criteria

    def setup_arcparam(self, arcparam, disperser=None, **null_kwargs):
        """
        Setup the arc parameters

        Args:
            arcparam: dict
            disperser: str, REQUIRED
            **null_kwargs:
              Captured and never used

        Returns:
            arcparam is modified in place

        """
        arcparam['lamps'] = ['CuI', 'ArI', 'ArII'] #  May be a handful of CuII lines too
        arcparam['nonlinear_counts'] = self.detector[0]['nonlinear']*self.detector[0]['saturation']
        if 'R150' in disperser:
            arcparam['n_first']=2 # Too much curvature for 1st order
            arcparam['disp']=0.63 # Ang per pixel (unbinned)
            arcparam['b1']= 4.54698031e-04
            arcparam['b2']= -6.86414978e-09
            arcparam['wvmnx'][1] = 6000.
            arcparam['wv_cen'] = 4000.
        elif 'R400' in disperser:
            arcparam['disp']=0.67 # Ang per pixel (unbinned) :: E2V  --  Hamamatsu is 0.74
            arcparam['min_ampl'] = 1000.0
        elif 'B600' in disperser:
            arcparam['n_first']=2 # Too much curvature for 1st order
            arcparam['disp']=0.63 # Ang per pixel (unbinned)
            arcparam['b1']= 4.54698031e-04
            arcparam['b2']= -6.86414978e-09
            arcparam['wvmnx'][0] = 3800.
            arcparam['wvmnx'][1] = 8000.
            arcparam['wv_cen'] = 4000.
        else:
            msgs.error('Not ready for this disperser {:s}!'.format(disperser))


class GeminiGMOSSSpectrograph(GeminiGMOSSpectrograph):
    """
    Child to handle Gemini/GMOS-S instrument
    """
    def __init__(self):

        # Get it started
        super(GeminiGMOSSSpectrograph, self).__init__()
        self.telescope = telescopes.GeminiSTelescopePar()
        self.spectrograph = 'gemini_gmos_south'
        self.camera = 'GMOS-S'
        self.detector = [
            # Detector 1
            DetectorPar(dataext         = 1,  # Not sure this is used
                        dispaxis        = 1,  # I think this is ignored, even if true
                        xgap            = 0.,
                        ygap            = 0.,
                        ysize           = 1.,
                        platescale      = 0.080,
                        darkcurr        = 0.0,
                        saturation      = 65535.,
                        nonlinear       = 0.86,
                        numamplifiers   = 4,
                        gain            = [1.83]*4,
                        ronoise         = [3.98]*4,
                        suffix          = '_01'
                        ),
            # Detector 2
            DetectorPar(dataext         = 2,  # Not sure this is used
                        dispaxis        = 1,
                        xgap            = 0.,
                        ygap            = 0.,
                        ysize           = 1.,
                        platescale      = 0.080,
                        darkcurr        = 0.0,
                        saturation      = 65535.,
                        nonlinear       = 0.86,
                        numamplifiers   = 4,
                        gain            = [1.83]*4,
                        ronoise         = [3.98]*4,
                        suffix          = '_02'
                        ),
            # Detector 3
            DetectorPar(dataext         = 3,  # Not sure this is used
                        dispaxis        = 1,
                        xgap            = 0.,
                        ygap            = 0.,
                        ysize           = 1.,
                        platescale      = 0.080,
                        darkcurr        = 0.0,
                        saturation      = 65535.,
                        nonlinear       = 0.86,
                        numamplifiers   = 4,
                        gain            = [1.83]*4,
                        ronoise         = [3.98]*4,
                        suffix          = '_03'
                        ),
        ]
        self.numhead = 13

    def header_keys(self):
        head_keys = self.gemini_header_keys()
        return head_keys

class GeminiGMOSNSpectrograph(GeminiGMOSSpectrograph):
    """
    Child to handle Gemini/GMOS-N instrument
    """
    def __init__(self):

        # Get it started
        super(GeminiGMOSNSpectrograph, self).__init__()
        self.telescope = telescopes.GeminiNTelescopePar()
        self.camera = 'GMOS-N'

    def header_keys(self):
        head_keys = self.gemini_header_keys()
        return head_keys

    def get_match_criteria(self):
        return self.gemini_get_match_criteria()


class GeminiGMOSNHamSpectrograph(GeminiGMOSNSpectrograph):
    """
    Child to handle Gemini/GMOS-N instrument with Hamamatsu detector
    """
    def __init__(self):

        # Get it started
        super(GeminiGMOSNHamSpectrograph, self).__init__()
        self.spectrograph = 'gemini_gmos_north_ham'

        self.detector = [  #  Hamamatsu (since 2011)
            # Detector 1
            DetectorPar(dataext         = 1,  # Not sure this is used
                        dispaxis        = 1,  # I think this is ignored, even if true
                        xgap            = 0.,
                        ygap            = 0.,
                        ysize           = 1.,
                        platescale      = 0.0807,
                        darkcurr        = 0.0,
                        saturation      = 129000.,
                        nonlinear       = 0.95,
                        numamplifiers   = 4,
                        gain            = [1.63]*4,
                        ronoise         = [4.14]*4,
                        suffix          = '_01'
                        ),
            # Detector 2
            DetectorPar(dataext         = 2,  # Not sure this is used
                        dispaxis        = 1,
                        xgap            = 0.,
                        ygap            = 0.,
                        ysize           = 1.,
                        platescale      = 0.0807,
                        darkcurr        = 0.0,
                        saturation      = 123000.,
                        nonlinear       = 0.95,
                        numamplifiers   = 4,
                        gain            = [1.63]*4,
                        ronoise         = [4.14]*4,
                        suffix          = '_02'
                        ),
            # Detector 3
            DetectorPar(dataext         = 3,  # Not sure this is used
                        dispaxis        = 1,
                        xgap            = 0.,
                        ygap            = 0.,
                        ysize           = 1.,
                        platescale      = 0.0807,
                        darkcurr        = 0.0,
                        saturation      = 125000.,
                        nonlinear       = 0.95,
                        numamplifiers   = 4,
                        gain            = [1.63]*4,
                        ronoise         = [4.14]*4,
                        suffix          = '_03'
                        ),
        ]
        self.numhead = 13


class GeminiGMOSNE2VSpectrograph(GeminiGMOSNSpectrograph):
    """
    Child to handle Gemini/GMOS-N instrument with E2V detector
    """
    def __init__(self):

        # Get it started
        super(GeminiGMOSNE2VSpectrograph, self).__init__()

        self.spectrograph = 'gemini_gmos_north_e2v'

        self.detector = [  #  E2V
            # Detector 1
            DetectorPar(dataext         = 1,  # Not sure this is used
                        dispaxis        = 0,  # I think this is ignored
                        xgap            = 0.,
                        ygap            = 0.,
                        ysize           = 1.,
                        platescale      = 0.0728,  # arcsec per pixel
                        darkcurr        = 0.0,
                        saturation      = 110900.,
                        nonlinear       = 0.95,
                        numamplifiers   = 2,
                        gain            = [2.27]*2,
                        ronoise         = [3.32]*2,
                        suffix          = '_01'
                        ),
            # Detector 2
            DetectorPar(dataext         = 2,  # Not sure this is used
                        dispaxis        = 0,
                        xgap            = 0.,
                        ygap            = 0.,
                        ysize           = 1.,
                        platescale      = 0.0728,
                        darkcurr        = 0.0,
                        saturation      = 115500.,
                        nonlinear       = 0.95,
                        numamplifiers   = 2,
                        gain            = [2.27]*2,
                        ronoise         = [3.32]*2,
                        suffix          = '_02'
                        ),
            # Detector 3
            DetectorPar(dataext         = 3,  # Not sure this is used
                        dispaxis        = 0,
                        xgap            = 0.,
                        ygap            = 0.,
                        ysize           = 1.,
                        platescale      = 0.0728,
                        darkcurr        = 0.0,
                        saturation      = 116700.,
                        nonlinear       = 0.95,
                        numamplifiers   = 2,
                        gain            = [2.27]*2,
                        ronoise         = [3.32]*2,
                        suffix          = '_03'
                        ),
        ]
        self.numhead = 7

    def header_keys(self):
        head_keys = self.gemini_header_keys()
        head_keys[0]['exptime'] = 'EXPOSURE'
        return head_keys

def read_gmos(raw_file, det=1):
    """
    Read the GMOS data file

    Parameters
    ----------
    raw_file : str
      Filename
    detector_par : ParSet
      Needed for numamplifiers if not other things
    det : int, optional
      Detector number; Default = 1

    Returns
    -------
    array : ndarray
      Combined image 
    header : FITS header
    sections : list
      List of datasec, oscansec, ampsec sections
    """

    # Check for file; allow for extra .gz, etc. suffix
    fil = glob.glob(raw_file+'*') 
    if len(fil) != 1:
        msgs.error("Found {:d} files matching {:s}".format(len(fil)))

    # Read
    msgs.info("Reading GMOS file: {:s}".format(fil[0]))
    hdu = fits.open(fil[0])
    head0 = hdu[0].header
    head1 = hdu[1].header

    # Number of amplifiers (could pull from DetectorPar but this avoids needing the spectrograph, e.g. view_fits)
    numamp = (len(hdu)-1)//3

    # Setup for datasec, oscansec
    dsec = []
    osec = []

    # get the x and y binning factors...
    binning = head1['CCDSUM']
    xbin, ybin = [int(ibin) for ibin in binning.split(' ')]

    # First read over the header info to determine the size of the output array...
    datasec = head1['DATASEC']
    x1, x2, y1, y2 = np.array(parse.load_sections(datasec, fmt_iraf=False)).flatten()
    biassec = head1['BIASSEC']
    b1, b2, b3, b4 = np.array(parse.load_sections(biassec, fmt_iraf=False)).flatten()
    nxb = b2-b1 + 1

    # determine the output array size...
    nx = (x2-x1+1)*numamp + nxb*numamp
    ny = y2-y1+1

    # allocate output array...
    array = np.zeros( (nx, ny) )

    if numamp == 2:
        if det == 1: # BLUEST DETECTOR
            order = range(6,4,-1)
        elif det == 2: # BLUEST DETECTOR
            order = range(3,5)
        elif det == 3: # BLUEST DETECTOR
            order = range(1,3)
    elif numamp == 4:
        if det == 1: # BLUEST DETECTOR
            order = range(12,8,-1)
        elif det == 2: # BLUEST DETECTOR
            order = range(8,4,-1)
        elif det == 3: # BLUEST DETECTOR
            order = range(4,0,-1)
    else:
        debugger.set_trace()

    # insert extensions into master image...
    for kk, jj in enumerate(order):

        # grab complete extension...
        data, overscan, datasec, biassec, x1, x2 = gemini_read_amp(hdu, jj)
                            #, linebias=linebias, nobias=nobias, $
                            #x1=x1, x2=x2, y1=y1, y2=y2, gaindata=gaindata)
        # insert components into output array...
        inx = data.shape[0]
        xs = inx*kk
        xe = xs + inx

        # insert data...
        # Data section
        #section = '[:,{:d}:{:d}]'.format(xs, xe)  # Eliminate lines
        section = '[{:d}:{:d},:]'.format(xs, xe)  # Eliminate lines
        dsec.append(section)
        array[xs:xe, :] = np.flipud(data)

        #; insert postdata...
        xs = nx - numamp*nxb + kk*nxb
        xe = xs + nxb
        #debugger.set_trace()
        #section = '[:,{:d}:{:d}]'.format(xs, xe)
        osection = '[{:d}:{:d},:]'.format(xs, xe)  # TRANSPOSED FOR WHAT COMES
        osec.append(osection)
        array[xs:xe, :] = overscan

    # make sure BZERO is a valid integer for IRAF
    obzero = head1['BZERO']
    #head0['O_BZERO'] = obzero
    head0['BZERO'] = 32768-obzero

    # Return, transposing array back to goofy Python indexing
    return array, head0, (dsec, osec)


def gemini_read_amp(inp, ext):
    """
    Read one amplifier of an LRIS multi-extension FITS image

    Parameters
    ----------
    inp: tuple 
      (str,int) filename, extension
      (hdu,int) FITS hdu, extension

    Returns
    -------
    data
    predata
    postdata
    x1
    y1

    ;------------------------------------------------------------------------
    function lris_read_amp, filename, ext, $
      linebias=linebias, nobias=nobias, $
      predata=predata, postdata=postdata, header=header, $
      x1=x1, x2=x2, y1=y1, y2=y2, GAINDATA=gaindata
    ;------------------------------------------------------------------------
    ; Read one amp from LRIS mHDU image
    ;------------------------------------------------------------------------
    """
    # Parse input
    if isinstance(inp, str):
        hdu = fits.open(inp)
    else:
        hdu = inp

    head1 = hdu[1].header

    # Deal with binning
    binning = head1['CCDSUM']
    xbin, ybin = [int(ibin) for ibin in binning.split(' ')]

    # get entire extension...
    temp = hdu[ext].data.transpose()
    tsize = temp.shape
    nxt = tsize[0]

    # parse the DETSEC keyword to determine the size of the array.
    header = hdu[ext].header
    detsec = header['DETSEC']
    x1, x2, y1, y2 = np.array(parse.load_sections(detsec, fmt_iraf=False)).flatten()

    # parse the DATASEC keyword to determine the size of the science region (unbinned)
    datasec = header['DATASEC']
    xdata1, xdata2, ydata1, ydata2 = np.array(parse.load_sections(datasec, fmt_iraf=False)).flatten()

    # grab the components...
    data = temp[xdata1-1:xdata2,:]

    # Overscan
    biassec = header['BIASSEC']
    xdata1, xdata2, ydata1, ydata2 = np.array(parse.load_sections(biassec, fmt_iraf=False)).flatten()
    overscan = temp[xdata1-1:xdata2,:]

    # Return
    return data, overscan, datasec, biassec, x1, x2



