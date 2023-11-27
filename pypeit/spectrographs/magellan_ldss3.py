"""
Module for Magellan/LDSS3 specific methods.

Important Notes:
    - 

"""

import numpy as np
import glob
from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.core import parse
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container
from pypeit import io

class MagellanLDSS3Spectrograph(spectrograph.Spectrograph):
    """
    Child to handle Magellan/LDSS3 specific code
    """

    name = 'magellan_ldss3'
    ndet = 1
    telescope = telescopes.MagellanTelescopePar()
    camera = 'LDSS3-C'
    url = 'https://www.lco.cl/technical-documentation/ldss-3-user-manual/'
    header_name = 'LDSS3-C'

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the ``PypeIt``-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """

        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card='RA')
        self.meta['dec'] = dict(ext=0, card='DEC')
        self.meta['target'] = dict(ext=0, card='OBJECT')
        self.meta['decker'] = dict(ext=0, card='APERTURE') # tells us about the slit used
        self.meta['dichroic'] = dict(ext=0, card=None, default='default')
        self.meta['binning'] = dict(card=None, compound=True) # ext=0, card='BINNING', default='1,1')

        self.meta['mjd'] = dict(card=None, compound=True) # ext=0, card='JD') # also have various times in the header, but this only jd
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=0, card='GRISM')
        self.meta['idname'] = dict(ext=0, card='EXPTYPE') 
        self.meta['instrument'] = dict(ext=0, card='INSTRUME')
        self.meta['amp'] = dict(ext=0, card='OPAMP') 

    def compound_meta(self, headarr, meta_key):
        """
        Methods to generate metadata requiring interpretation of the header
        data, instead of simply reading the value of a header card.
        
        Args:
            headarr (:obj:`list`):
                List of `astropy.io.fits.Header`_ objects.
            meta_key (:obj:`str`):
                Metadata keyword to construct.

        Returns:
            object: Metadata value read from the header(s).
        """

        if meta_key == 'binning':
            
            binspatial, binspec = parse.parse_binning(headarr[0]['BINNING'])
            binning = parse.binning2string(binspec, binspatial)
            return binning

        elif meta_key == 'mjd':

            jd = headarr[0]['JD']
            mjd = jd - 2400000.5
          
            return mjd

    def pypeit_file_keys(self):
        """
        Define the list of keys to be output into a standard ``PypeIt`` file.

        Returns:
            :obj:`list`: The list of keywords in the relevant
            :class:`~pypeit.metadata.PypeItMetaData` instance to print to the
            :ref:`pypeit_file`.
        """
        pypeit_keys = super().pypeit_file_keys()
        # TODO: Why are these added here? See
        # pypeit.metadata.PypeItMetaData.set_pypeit_cols
        pypeit_keys += ['calib', 'comb_id', 'bkg_id']
        return pypeit_keys


    def get_detector_par(self, hdu, det):
        """
        Return metadata for the selected detector.

        Args:
            det (:obj:`int`):
                1-indexed detector number.
            hdu (`astropy.io.fits.HDUList`_, optional):
                The open fits file with the raw image of interest.  If not
                provided, frame-dependent parameters are set to a default.

        Returns:
            :class:`~pypeit.images.detector_container.DetectorContainer`:
            Object with the detector metadata.
        """

        # get binning from the meta value (which is compound, i.e. parsed from the header, but needs editing)
        binning = self.get_meta_value(self.get_headarr(hdu), 'binning')

        # Detector 1
        detector_dict = dict(
            binning         = binning,
            det             = 1,
            dataext         = 1,
            specaxis        = 1,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.189,
            darkcurr        = 25.0, # e-/pix/hr -- are these the correct units?
            saturation      = 65536., # ADU
            nonlinear       = 1.0, # have not touches this
            mincounts       = -1e10,
            numamplifiers   = 1, # from one of the headers
            gain            = np.atleast_1d(1.67), # fast mode from user manual for C1
            ronoise         = np.atleast_1d(7.0),
            datasec         = np.atleast_1d('[1:1024,1:4096]'),
            oscansec        = np.atleast_1d('[1025:1152,1:4096]')
            )
        return detector_container.DetectorContainer(**detector_dict)

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.
        
        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of ``PypeIt`` methods.
        """
        par = super().default_pypeit_par()

        # Wavelengths
        # 1D wavelength solution
        par['calibrations']['wavelengths']['rms_threshold'] = 0.5
        par['calibrations']['wavelengths']['sigdetect'] = 5
        par['calibrations']['wavelengths']['lamps'] = ['ArI','NeI','HeI']
        par['calibrations']['wavelengths']['reid_arxiv'] = 'wvarxiv_magellan_ldss3_20231122T0943.fits' 
        par['calibrations']['wavelengths']['method'] = 'full_template' # 'holy-grail' # if things fail with this, try reidentify (but note that archived wavelength solution is required)

        # Set slits and tilts parameters
        par['calibrations']['tilts']['tracethresh'] = 20. # default value
        par['calibrations']['slitedges']['sync_predict'] = 'nearest' # required for long-slit

        # Do not require dark frames or overscans
        turn_off = dict(use_darkimage=False, use_overscan=False)
        par.reset_all_processimages_par(**turn_off)

        # LACosmics parameters
        par['scienceframe']['process']['mask_cr'] = True
        par['scienceframe']['process']['sigclip'] = 4.0
        par['scienceframe']['process']['objlim'] = 1.5

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['standardframe']['exprng'] = [None, 900]
        par['calibrations']['arcframe']['exprng'] = [0, 12]
        par['calibrations']['darkframe']['exprng'] = [0, None]
        par['scienceframe']['exprng'] = [20, None]

        # Sensitivity function parameters
        par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['polyorder'] = 5
        par['sensfunc']['IR']['telgridfile'] = 'TelFit_LasCampanas_3100_26100_R20000.fits'

        return par

    def bpm(self, filename, det, shape=None, msbias=None):
        """
        Generate a default bad-pixel mask.

        Even though they are both optional, either the precise shape for
        the image (``shape``) or an example file that can be read to get
        the shape (``filename`` using :func:`get_image_shape`) *must* be
        provided.
        Args:
            filename (:obj:`str` or None):
                An example file to use to get the image shape.
            det (:obj:`int`):
                1-indexed detector number to use when getting the image
                shape from the example file.
            shape (tuple, optional):
                Processed image shape
                Required if filename is None
                Ignored if filename is not None
            msbias (`numpy.ndarray`_, optional):
                Master bias frame used to identify bad pixels
        Returns:
            `numpy.ndarray`_: An integer array with a masked value set
            to 1 and an unmasked value set to 0.  All values are set to
            0.
        """
        # Get the empty bpm: force is always True
        bpm_img = self.empty_bpm(filename, det, shape=shape)

        # Fill in bad pixels if a master bias frame is provided
        if msbias is not None:
            return self.bpm_frombias(msbias, det, bpm_img)

        return bpm_img 

    def configuration_keys(self):
        """
        Return the metadata keys that define a unique instrument
        configuration.

        This list is used by :class:`~pypeit.metadata.PypeItMetaData` to
        identify the unique configurations among the list of frames read
        for a given reduction.

        Returns:
            :obj:`list`: List of keywords of data pulled from file headers
            and used to constuct the :class:`~pypeit.metadata.PypeItMetaData`
            object.
        """
        return ['dispname', 'binning'] # not really sure here what's right for LDSS3

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.

        Args:
            ftype (:obj:`str`):
                Type of frame to check. Must be a valid frame type; see
                frame-type :ref:`frame_type_defs`.
            fitstbl (`astropy.table.Table`_):
                The table with the metadata for one or more frames to check.
            exprng (:obj:`list`, optional):
                Range in the allowed exposure time for a frame of type
                ``ftype``. See
                :func:`pypeit.core.framematch.check_frame_exptime`.

        Returns:
            `numpy.ndarray`_: Boolean array with the flags selecting the
            exposures in ``fitstbl`` that are ``ftype`` type frames.
        """
        # science, arcs are "Object", bias are "Bias", flats are "Flat"
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            return good_exp & (fitstbl['idname'] == 'Flat')
        elif ftype == 'bias':
             return good_exp & (fitstbl['idname'] == 'Bias')
        elif ftype == 'standard':
            return good_exp & (fitstbl['idname'] == 'Object')
        elif ftype == 'science':
            return good_exp & (fitstbl['idname'] == 'Object')
        elif ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['idname'] == 'Object')
        else:
            msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
            return np.zeros(len(fitstbl), dtype=bool)

    def get_rawimage(self, raw_file, det):
        """
        Read raw images and generate a few other bits and pieces
        that are key for image processing.

        Parameters
        ----------
        raw_file : :obj:`str`
            File to read
        det : :obj:`int`
            1-indexed detector to read

        Returns
        -------
        detector_par : :class:`pypeit.images.detector_container.DetectorContainer`
            Detector metadata parameters.
        raw_img : `numpy.ndarray`_
            Raw image for this detector.
        hdu : `astropy.io.fits.HDUList`_
            Opened fits file
        exptime : :obj:`float`
            Exposure time read from the file header
        rawdatasec_img : `numpy.ndarray`_
            Data (Science) section of the detector as provided by setting the
            (1-indexed) number of the amplifier used to read each detector
            pixel. Pixels unassociated with any amplifier are set to 0.
        oscansec_img : `numpy.ndarray`_
            Overscan section of the detector as provided by setting the
            (1-indexed) number of the amplifier used to read each detector
            pixel. Pixels unassociated with any amplifier are set to 0.
        """

#        # Open
#        hdu = io.fits_open(raw_file)
#
#        # Grab the DetectorContainer and extract the raw image
#        detector = self.get_detector_par(det=det, hdu=hdu)
#        raw_img = hdu[detector['dataext']].data.astype(float)
#
#        # Exposure time (used by RawImage) from the header
#        headarr = self.get_headarr(hdu)
#        exptime = self.get_meta_value(headarr, 'EXPTIME')
#
#        for section in ['DATASEC', 'BIASSEC']:
#            # Get the data section from Detector
#            image_sections = detector[section]
#
#            # Initialize the image (0 means no amplifier)
#            pix_img = np.zeros(raw_img.shape, dtype=int)
#            for i in range(detector['numamplifiers']):
#
#                if image_sections is not None:
#                    # Convert the (FITS) data section from a string to a slice
#                    # DO NOT send the binning (default: None)
#                    datasec = parse.sec2slice(image_sections[i], one_indexed=True,
#                                              include_end=True, require_dim=2)
#                    # Assign the amplifier
#                    pix_img[datasec] = i+1
#
#            # Finish
#            if section == 'DATASEC':
#                rawdatasec_img = pix_img.copy()
#            else:
#                oscansec_img = pix_img.copy()
#
#        # Return
#        return detector, raw_img, hdu, exptime, rawdatasec_img, oscansec_img

        # Check for file; allow for extra .gz, etc. suffix
        fil = glob.glob(raw_file + '*')
        if len(fil) != 1:
            msgs.error("Found {:d} files matching {:s}".format(len(fil), raw_file))

        # load in the data
        hdu = io.fits_open(fil[0])
        head1 = hdu[0].header

        # read header info to determine size
        datasec = head1['DATASEC']
        x1, x2, y1, y2 = np.array(parse.load_sections(datasec, fmt_iraf=False)).flatten()
        biassec = head1['BIASSEC']
        b1, b2, b3, b4 = np.array(parse.load_sections(biassec, fmt_iraf=False)).flatten()
        nxb = b2 - b1 + 1 # number of x pixels in the bias scan
        nx = (x2 - x1 + 1) + nxb # number of x pixels in the array
        ny = y2 - y1 + 1 # number of y pixels in the array

        # create arrays to save
        array = hdu[0].data.T[:, :ny] * 1.0
        data = array[:nx-nxb,:] # remove the bias scan region from data
        overscan = array[nx-nxb:,:] # overscan region is just the bias region

        # get detector parameters with method we wrote above
        detector_par = self.get_detector_par(hdu, det)

        # for if we were going to use both arrays
        #data1, overscan1, datasec1, biassec1, x1_1, x2_1, nxb1 = ldss3_read_amp(fil[0])
        #data2, overscan2, datasec2, biassec2, x1_2, x2_2, nxb2 = ldss3_read_amp(fil2[0])
        #nx, ny = x2_1 + x2_2 + nxb1 + nxb2, data1.shape[1]

        # allocate output array...
        array = np.zeros((nx, ny))
        rawdatasec_img = np.zeros_like(array, dtype=int)
        oscansec_img = np.zeros_like(array, dtype=int)

        ## For amplifier 1
        array[:nxb,:] = overscan # set the overscan region to have the overscan data in it
        array[nxb:,:] = data # set the data region to have the data in it
        rawdatasec_img[nxb:, :] = 1 # set where the data is to 1
        oscansec_img[:nxb,:] = 1 # set where the overscan region is to 1

        ## transpose so that it picks out the right edges
        #array = array.T
        #rawdatasec_img = rawdatasec_img.T
        #oscansec_img = oscansec_img.T

        # Need the exposure time
        exptime = hdu[self.meta['exptime']['ext']].header[self.meta['exptime']['card']]

        return detector_par, array, hdu, exptime, rawdatasec_img, oscansec_img
