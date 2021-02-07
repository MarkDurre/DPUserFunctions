# QFitsView DPUser
## Setting Up DPUSER Functions Library
1. Place all libraries (files with “lib\_\*.dpuser” name) in a convenient location  e.g. “*my\_code\_location*/DPUserlib/Functions" (subsituting for *my_code_location* e.g. "/Users/*user*/Programs").

2. Place “startup.dpuser” in e.g. “*my\_code\_location*/DPUserlib". This assumes that all the library functions are located in the "Function" sub-folder and runs the "lib\_all.dpuser" script. This file must be modified for your own requirements; it also sets the *DPUSER_DIR* environment variable that can be accessed by other scripts, using the **getenv** function in QFitsView.

    ```
    //Example startup.dpuser
    //This line must be modified (change *my_code_location*) for individual users
    setenv "DPUSER_DIR", "*my_code_location*/dpuserlib"
    //dpuserdir is used by libraries to call other libraries if required
    dpuserdir=getenv("DPUSER_DIR")
    print "Running General Functions : DPUser Directory - "+dpuserdir
    run dpuserdir+"/functions/lib_all.dpuser"
    print "Finished General Functions - "+dpuserdir
    //You can put your own startup dpuser code here
    ```

3. Create a symbolic link in the root directory to this folder (for macOS 10.15+ use the `/etc/synthetic.conf` symbolic links method - reference [here](#https://stackoverflow.com/questions/58396821/what-is-the-proper-way-to-create-a-root-sym-link-in-catalina)). This link must be called "dpuserlib".

4. When QFitsView starts, it automatically looks for and executes the script file "/dpuserlib/startup.dpuser" i.e. it runs the script set up above. This runs all libraries to make the functions available to QFitsView - you will see a whole bunch of “Stored function…” and “Stored procedure…” plus “Finished General Functions” text lines on the DPUSER area in QFitsView.

5. If the above folder conventions are not used, the following files must be modified:
    - *my_program_location*/DPUserlib/startup.dpuser
    - *my_program_location*/DPUserlib/Functions/lib\_all.dpuser

##  Application Global Variables
These are defined internally; they can be overwritten in a QFitsView session.

c          =     299792458.0<br>
pi         =     3.14159<br>
e          =     2.71828<br>
naxis1     =     256<br>
naxis2     =     256<br>
plotdevice =     /QT<br>
method     =     0<br>
tmpmem     =     20971520<br>

## Editing DPUser Code

The main documentation for DPUser is through [this link](#https://www.mpe.mpg.de/~ott/dpuser/). DPUser code can be edited with QFitsView *DPUSER > Script Editor*. It can also be edited by various external text editors; e.g. **BBEdit**. To facilitate this, a language module has been implemented - "DPUser.plist", with the following highlighting features: 

- Syntax - both structural commands - e.g. "if", "else" etc. and internal DPUser functions/procedures. As new functions/procdures are implemented in QFitsView, the "BBLMPredefinedNameList" array must be updated.
- Comments (both for "/\*..\*/" and "//").
- Strings ("..")
- Function and procedure prefixes

The file is installed in BBEdit's language module directory - by default "/Users/*username*/Library/Application Support/BBEdit/Language Modules".

Note that after editing code in an external text editor, the script must be executed again with QFitsView.

##  Parameter data types
* FITS or buffer input - *cube* is 3D, *image* is 2D, *spectrum* is 1D, *data* is 1, 2 or 3D, *mef* is multi-extension FITS data (created using the **list** function for up to 10 data buffers).
	* Pixel co-ordinates
	* p, p1, p2,… (general pixels co-ordinates)
	* x, x1, x2,…, y, y1, y2,…, z, z1, z2,… (for x, y or z axes)
	* pixel/wavelength masking - an array of 2n values [x1,x2,x3,x4...] - paired pixel numbers or wavelengths (pairs x1..x2,x3..x4 etc.) for spectral masking
* WCS co-ordinates
   * wcs - for axis set [CRPIXn, CRVALn, CDELTn]
  * w, w1, w2,…. (for individual coordinates)
  * l, l1, l2, …. (wavelengths)
* Maps
   * velmap - velocity map cube (either standard *velmapstd* or extended *velmapext* format)
  * wvtmap - weighted Voronoi tessellation map (region numbers)
## Libraries
### lib_all

Runs all the following libraries; this just consists of script lines to execute other scripts in the "Functions" sub-folder. Other scripts can be executed by adding appropriate lines, e.g.

`@SomeFolder/SomeScript.dpuser`

If full path to script is not given, it is assumed to be relative to the DPUser folder defined above.

### lib_wcs

**<u>*Transform to and from World Coordinate Systems and Pixels*</u>**

**function get_WCS_data, data, axis** - return array [CRVAL, CRPIX, CDELT] for *axis* (1,2 or 3) of *data*

**procedure set_WCS_data, data, wcs, axis** - sets WCS values for *data* for *axis* (1,2 or 3)

**function cvt_pixel_WCS, pix, crpix, crval, cdelt** - convert pixel number *pix* to WCS coordinates using *crpix, crval, cdelt*, asssuming linear projection (cartesian). The DPUser **worldpos** function does this including the correct projections.

**function cvt_WCS_pixel, value, crpix, crval, cdelt** - convert WCS coordinate *value* to pixel using *crpix, crval, cdelt*, asssuming linear projection (cartesian).  The DPUser **pixpos** function does this including the correct projections.

**function cvt_WCS_pixel_data, data, value, axis** - converts pixel *value* to WCS for *data* for *axis*

**function cvt_pixel_WCS_data, data, value, axis** - reverse of *cvt_WCS_pixel_data*

**function WCS_range, data, p1, p2, axis, prntflag** - returns WCS coordinates as range [*w1,w2*] from pixel values [*p1,p2*] for axis on *data*, if *prntflag*=1, then print range

**function pixel_range, data, w1, w2, axis, prntflag**  - returns pixel values as range [*p1,p2*] from WCS co-ordinates [*w1,w2*] for axis on *data*, if *prntflag*=1, then print range

**function set_WCS_default, data** - checks *data* has minimal WCS keys set (to 1 by default).

<a name="get_WCS_image"></a>**function get_WCS_image, image** - Get axis 1 and 2 WCS data for *image* (can be cube) and calculates CDELT and rotation angle from CD keys; result is array of key values [CRPIX1, CRVAL1, CD1_1, CD1_2, CRPIX2, CRVAL2, CD2_1, CD2_2, CDELT1, CDELT2, CROTA2]

**function set_WCS_image_scale, image, wcs2d, xscale, yscale** - Rescales (e.g. for non-integer re-binning) using *xscale* and *yscale* and sets WCS data for *image* (or cube), including CD keys - removes CDELT and CROTA2 keywords. Input *wcs2d* is format as for **[get_WCS_image](#get_WCS_image)**.

**function get_fits_key, data, key** - substitute for **getfitskey** with check that *key* exists; returns blank if not found.

**function get_WCS_values, data** - create WCS array [1,cv,cd] from dispersion *data*, calculated from first/last values and number of elements.

**function WCS_cdelt_cd, cdelt1, cdelt2, rotang** - converts WCS x,y pixel sizes (*cdelt1, cdelt2*) and rotation angle (*rotang*) to CD matrix values. Returns vector [CD11, CD12, CD21, CD22].

**procedure set_cd_keys, data, cdkeys** - sets CD keys in FITS header of *data* from vector *cdkeys* (in same format as **WCS_cdelt_cd** function) and deletes the CDELT1/2 keys.

### lib_cube
**<u>*Data cube functions*</u>**

**function cube_trim_xy, cube, x1, x2, y1, y2**- sets cube to zero for x<*x1*, x>*x2*, y<*y1*, y>*y2* (**cblank** cube first)

**function cube_trim_wl, cube, w1, w2, value, trimflag** - sets *cube* to *value* (usually zero) for l < *l1*, l > *l2* in axis 3 using WCS (cblank cube first). If *trimflag*=1, then truncate the cube outside the wavelength range.

**function cube_spectrum_mask, cube, mask, level** - mask *cube* on spectral wavelength with *mask* (pixel pairs), set masked pixels to *level*

**function cube_clip, cube, lvl, thresh, mask** - clips *cube* <0 and > *lvl* in image (x-y) plane, does dpixapply using threshold, cube is spectrally masked

**function cube_clip_y, cube, lvl, thresh** - as above, but in the x-z plane    

**function cube_interp_z, cube, x1, x2, y1, y2, z1, z2** - interpolate *cube* in image plane over rectangle [*x1:x2,y1:y2*] in each of wavelength plane [*z1:z2*]

<a name="cube_interp_x"></a>**function cube_interp_x, cube, x1, x2, y1, y2, z1,z2** - interpolate *cube* in image plane over rectangle [*y1:y2,z1:z2*] in each of spatial range [*x1:x2*]

<a name="cube_interp_y"></a>**function cube_interp_y, cube, x1, x2, y1, y2, z1,z2** - interpolate *cube* in image plane over rectangle [*x1:x2,z1:z2*] in each of spatial range [*y1:y2*]

<a name="cube_interp_xy"></a>**function cube_interp_xy, cube, x1, x2, y1, y2, z1,z2** - as above, but interpolate over wavelength [*z1:z2*] in the xy plane

**function cube_set_value, cube, x1, x2, y1, y2, z, xv, yv** - set rectangle [*x1:x2,y1:y2*] at image plane *z* to value at [*xv, yv*]

**function cube_pixfix_xy, cube, pixfixdata, n** - fix cube using **cube_interp_xy** function, n sets, *pixfixdata* is  *n* x 6 array [x1,x2,y1,y2,z1,z2].

**function cube_single_pixel_fix, cube, x, y** - do **cube_interp_xy** for all z axis for single spaxel

**procedure cube_bit_nan, cube,x,y**  - set spaxel [*x*,*y*] to 0/0 along whole *cube* z axis

**function cube_clean_dpix, cube, scale** - Clean *cube*  by dpixcreate/apply, the threshold for dpixcreate is set from maximum of median image divided by *scale*.

**function cube_resize_center, cube, xcen, ycen, xsize, ysize** - resize *cube*  to size *xsize*,*ysize* and center on pixel [*xcen,ycen*]

**function cube_shift_xy, cube, xshift, yshift** - sub-pixel shifts *cube* by [*xshift, yshift*]

**function cube_redisp, cube, disp_old, disp_new, prnt** - redisperses *cube* (axis 3) to new from *disp_old* to *disp_new* dispersion spectra by interpolation. If *prnt* <> 0, print diagnostics.

**function cube_symm_flip, cube, lambda, width, part** - symmetrically flip *cube* about wavelength *lambda*, *part* =0 (left) or 1 (right), trims cube to *lambda*+-*width*

**function cube_rotate, cube, xcen, ycen, rot_angle, pixscale** - rotate *cube* on center [*xcen*,*ycen*] by *rot_angle*, setting *pixscale* in arcsec/pixel. This is used because the **rotate** function does not work on cubes.

**function cube_centroids, cube** - get centroids at each wavelength layer (z axis). Returns a FITS array of dimensions naxis3(*cube*) x 2, with x and y centroids at each pixel layer. The wavelength WCS is set.

<a name="cube_centroids_gauss"></a>**function cube_centroids_gauss, cube, xe, ye, we, mask** - get centroid at each pixel layer, with estimated center at [*xe,ye*] over fitting window *we*. The spectrum is masked by *mask* (if not zero or not entered).

**function cube_centroid_gauss_align, cube, xc, yc, xe, ye, we, mask** - align *cube* centroids at each pixel layer, using [*xe, ye, we, mask*] are the estimate parameters of the peak (as for **[cube_centroids_gauss](#cube_centroids_gauss)**), with the centroids aligned to [*xc, yc*].

**function cube_cont_slope, cube, mask** - returns image with continuum slope at each spaxel of *cube*, masked by wavelength pairs *mask*

**function cube_spectrum_add, cube, spectrum, x1, x2, y1, y2** - adds *spectrum* to *cube* for each spaxel. If any of *x1, x2, y1, y2* are set to zero, then perform the action over the pixel range. By default (not given), these are set to zero. To subtract *spectrum*, add by negative.

**function cube_spectrum_multiply, cube, spectrum, x1, x2, y1, y2** - multiply *cube* by *spectrum* as for **cube_spectrum_add**. To divide by *spectrum*, multiply by inverse.

**function cube_set_pixlayers, cube, pixl, p1, p2** - set *cube* layers [*p1, p2*] to the values for layer *pixl*.

**function cube_wavelength_correct, cube, correction** - corrects the wavelength solution at each spaxel by shifting the spectrum.  *correction* is an image of the same dimensions as the *cube* x and y axes and is in same units as the spectral axis of *cube*. 

<a name="cube_velocity_correct"></a>**function cube_velocity_correct, cube, velmodel** - applies a velocity field model *velmodel* (in km/s) to each spaxel of *cube*. Each spaxel is re-dispersed by intepolation.

**function cube_to_2d, cube** - Convert data *cube* to 2d apertures for IRAF. Returns a 2D array with spectrum on the x-axis and all spaxels on the y axis.

**function cube_set_flags_nan, cube, layer** - set up flags image for cube_interp_flags, from a data *cube* (e.g. a velmap) from *layer*. This retuens an image with same dimensions as x and y axes as the cube, with 1 where pixel in “NaN”, 0 else.

**function cube_interp_flags, cube, flags, xi1, xi2, yi1, yi2, dmax** - interpolate over pixels in *cube* where *flags* is set to 1, 0 = good values to use for interpolation. [*xi1:xi2, yi1:yi2*] is region to interpolate (*xi1* = 0 - do whole area). *dmax* is maximum distance from “good” pixels. *flags* can be generated from **cube_set_flags_nan**.

**function cube_deslope, cube, mask, wlflag** - deslope *cube* for each spectrum using **spectrum_deslope**. *wlflag* = 1 if mask values in wavelength

**function cube_clean_pixels, cube, layer, npix** - Remove singleton pixels surrounded by Nan’s, opposite of **cube_interp_flags**, used to clean up boundaries etc. *npix* is max number of good pixels around each pixel before blanking.

**function cube_radial_spectrum, cube, xc, yc, rstep, nstep, ann** - Radial spectra of *cube*, centered [*xc,yc*] radial steps *rstep*, number of steps *nstep*. If *ann=1*, output annular spectra

**function cube_from_image_spectrum, image, spectrum** - Creates a cube from an image and spectrum. Wavelength axis of cube is spectrum scaled by image value.

**function cube_rebinxy, cube, xscale, yscale, kernel** - Rebin *cube* or image pixel scaling in x and y directions by *xscale*, *yscale*. Uses the **interpolate** dpuser function with kernel *kernel*. Note this function DOES NOT handle the WCS co-ordinates scaling; use **get_WCS_cube** and **set_WCS_cube_scale** functions.

**function cube_rebinfrac, inbuff, xscale, yscale** - Rebins *cube* (image) to *xscale*, *yscale* using fractional binning. Note comments about WCS values as above.

**function cube_rebinx, cube, xscale**  - Rebin *cube* (image) pixel scaling in x direction ONLY. Uses the **interpol** dpuser function (quicker than **interpolate**). Note comments about WCS values as above.

**function cube_combine_avg, mef, omit** - combine cubes by averaging. *mef* is a multi-extension fits list, *omit* is a value to reject in the averaging. If not provided or set to 0/0, then no rejection is done. *mef* is created from mutiple cubes (up to 10) by e.g. 
`mef = list(buffer1, buffer2, buffer3, ....)`

**function cube_combine_median, mef, omit** - median combine cubes as for **cube_combine_avg**.

<a name="cube_apply_snr"></a>**function cube_apply_snr, signal, snr, snrscale** - adds noise to a *signal* cube using signal-to-noise ratio *snr*, multiplied by factor *snrscale* (default 1). *snr* can be a single value, a spectrum which is applied at each spaxel, an image where a single value is applied at each spaxel or a cube with different values at each pixel. The noise is a random gaussian value with standard deviation of the applicable SNR.

### lib_image
**<u>*Image functions*</u>**

**function image_erodenan, image** - erode *image*, pixels set to Nan if any neighbour is Nan.

**function image_smooth, image, smooth** - smooth *image* with NaN values - *smooth* integer=boxcar, non-integer=gaussian.

**function image_interp_x, image, x1, x2, y1, y2** - as for **[cube_interp_xy](#cube_interp_xy)**, but for single *image* .

**function image_interp_y, image, x1, x2, y1, y2** - as for **[cube_interp_x](#cube_interp_x)**, but for single *image* .

**function image_interp_xy, image, x1, x2, y1, y2** - as for **[cube_interp_y](#cube_interp_y)**, but for single *image* .

**function image_from_profile, profile, xp, yp, xc, yc** - create 2D image from 1D *profile*, size of output image is *xp* x *yp* , [*xc, yc*] - center of rebuilt profile

**function image_bfilter, image, order, cutoff** - Butterworth filter an *image*, assume square image, filter order *=order*, *cutoff* =Nyquist cutoff (0-1)

**function image_enclosed_flux, inbuff, xc, yc, r, smth** - Get enclosed flux within radius *r* from [*xc, yc*] (pixels). If *smth*>0, Gaussian smooth the output 

**function image_avg, image, x, y, s** - average value of image in square aperture [*x,y*] +-*s* pixels.

**function image_structure, image, psf** - Returns structure map from *image* and *psf* by formula *image*/(*image* ⊗ *psf*) x *psf*^T, "⊗”=convolution, "^T" = transpose.

**function image_interp_flags, image, flags, xi1, xi2, yi1, yi2, dmax** - Interpolate *image* over flagged spaxels, *flags* - 2D data with same x/y axes size as image, with value=1 to be interpolated, value=0 - good pixels, [*xi1:xi2, yi1:yi2*] - co-ordinate range to interpolate over. If not input, then do all spaxels. *dmax* - maximum pixel distance for interpolation (=0 don't test)

**function image_cut, image, x, y, a** - does **twodcut** at [*x,y*] angle *a* and reset WCS correctly

### lib_spectrum
**<u>*Spectrum functions*</u>**

**function spectrum_make_disp, val, delt, pix, n** - make 1D vector over range defined by WCS *val*, *delt*, *pix*, *n*.

**function spectrum_make_disp_data, data, axis** - make 1D vector over range defined by WCS values from *data* axis (1,2 or 3)

**function spectrum_make_disp_n, val1, val2, n** - make 1D vector over range [*val1:val2*], number of points *n*

**function spectrum_mask, spectrum, mask, value, wlflag** - *spectrum* set to *value* between pixel pairs in *mask*. Works for 1D or 3D, assuming last axis is spectrum. *wlflag* =0, mask is in pixels, =1, mask is in wavelength. *value* is usually 0/0 to blank masked sections for polynomial fitting.

**function spectrum_cont_slope, spectrum, mask, wlflag** - continuum slope of *spectrum*, masked by wavelength *mask*/*wlflag* (as for **spectrum_mask**)

**function spectrum_deslope, spectrum, mask, wlflag** - deslope spectrum, using **spectrum_cont_slope** and *mask/wlflag* parameters

**function spectrum_polyfit, spectrum, order, mask, wlflag** - fit polynomial of *order* to  masked *spectrum* with *mask*, *wlflag*. Returns n x 3 array, 1st row=original data masked, 2nd row=polynomial fit, 3rd row = residual

**function spectrum_symm_flip, spectrum, lambda, part** - split *spectrum* at wavelength *lambda*, flip and add, taking left (*part*=0) or right (*part*=1) sections

**function spectrum_wave_to_lambda, spectrum, l1, l2, nl ** - converts a wavenumber *spectrum* to a wavelength spectrum. *l1*..*l2* are a wavelength range to interpolate over with *nl* points. By default, the wavelength range and number of points of the original spectrum are used. WCS values are set.

**function spectrum_wave_to_lambda, wndata** - convert wavenumber spectrum *wndata* to wavelength (nm) with same axis length

**function spectrum_make_gauss, spectrum, bi, bs, h, l, w** - make spectrum with gaussian from *spectrum* WCS. *bi, bs* - base intercept and slope, *h* -  height, *lc* - center wavelength, *w* - FWHM (creates artificial gaussian emission line).

**function spectrum_make_lorentz, spectrum, bi, bs, h, lc, w** - as for **spectrum_make_gauss** but makes a Lorentzian emission line.

**function spectrum_redisp_lin, spectrum, data, daxis, xmin,delt, npix, zero, norms,  prnt, fluxcons** - re-disperse a *spectrum*. Parameters are:

- *data* - data with dispersal solution (if =0 then use parameters for dispersion)
- *daxis* - spectral axis of data (default is last axis of *data*)
- *xmin, delt, npix* - dispersion solution if data=0
- *zero* - if =1, then set redispersed spectra to zero where out of original range, rather than NaN (=0)
- *norms* - if =1, normalize dispersed spectra [0,1] (default 0)
- *prnt*  - if =1, print spectral range information (default 0)
- *fluxcons*  - if = 1, spectrum is flux, rather than flux density, so conserve total (default 0)

**function spectrum_from_xy, spectrum** - re-disperse *spectrum* from 2D x and y bintable to wavelength range and same number of points.

**function spectrum_from_tablexy, data, l1, l2, npix, xscl, yscl** - re-disperse *spectrum* from 2D x and y bintable to wavelength range *l1..l2* and number of points *npix*. x and y values are scaled by *xscl* and *yscl* respectively.

**function spectrum_redisp, spectrum, l1, l2, npix, xscl, yscl** - as for **spectrum_from_tablexy** but from standard *spectrum*.

**function spectrum_from_dataxy, xdata, ydata, l1, l2, npix, xscl, yscl** - re-disperse spectrum from 2D x and y data to wavelength range [*w1, w2*] with step *delt*.

**function spectrum_interp, spectrum, x1, x2** - Smooth over bad pixels [*x1:x2*].

**function spectrum_sn, spectrum, window** - Estimate spectrum S/N from itself - not 100% accurate but good for comparisons, *window* is smoothing and noise estimation window. Returns vector of same length as *spectrum* with S/N estimate, blank where *spectrum* is 0.

**function spectrum_clean, inbuff, thresh** - clean *spectrum* using **dpixcreate**/**dpixapply**. If *thresh* is not sepecifed the threshold for **dpixcreate** is set to median(*spectrum*)/2.

**function spectrum_apply_snr, signal, snr** - as for **[cube_apply_snr](#cube_apply_snr)**. *snr* can be a single value or spectrum.

### lib_io
**<u>*Input/output to and from text and fits files*</u>**

**function io_text_FITS_1D, bintable2d** - converts string array buffer *bintable2d* with format of "wavelength, data" to spectrum fits data, setting WCS values. Assumes wavelength is evenly spaced. Note **import** function of QFitsView does very similar (with more parameters).

**function io_text_FITS_3D, bintable2d, nx, ny, nz, blank**  - converts string array buffer *bintable2d*, with format of “i,j,v1,v2..." to fits data cube size [*nx,ny,nz*]. Default value for resulting cube is *blank* (e.g. 0 or 0/0) - can have missing [*i,j*].

**function io_text_FITS_interp, fname, xstart, xdelta, xnum, xscale, yscale, ignore** - converts text from file *fname*, with format of "wavelength, data" to spectrum fits data, setting WCS values. The values are interpolated to the range defined by *xstart*, *xdelta* and *xnum*. Wavelength and data value are scaled by *xscale*, *yscale* (default 1). *Ignore* lines at the start are skipped (e.g. column headers).

<a name="io_FITS_text_1D"></a>**function io_FITS_text_1D, spectrum, prefix, cutoff** - converts *spectrum* to text, CSV format, line 1 = “*prefix*\__Wavelength, *prefix*\_Counts”. Values below *cutoff* (non-zero) are set to “NaN” (Be aware of QFitsView Edit > Copy functionality)

**function io_FITS_text_2D, image, prefix** - converts *image* to text, CSV format, line 1 = “*prefix*\_Wavelength, *prefix*\_Flux_1, *prefix*_Flux_2 …. "

**procedure io_FITS2TXT_1D, fname, cutoff** - converts 1D FITS to text file *fname* assuming file is in  working directory - output is same as input file with “.txt” type. *Cutoff* as for **[io_FITS_text_1D](#io_FITS_text_1D)**.

**procedure io_FITS2TXT_2D, fname** - as above but for image (2D) file

**function io_cube_from_xyz, cube,bintable2d, n** - make a cube from *bintable2d*, *cube* is template, resized to *n* on axis 3,  first 2 values in data are x,y co-ords, rest are values along z axis

**function io_import_TXT_1D, fname** - import data from file *fname* in text format

### lib_masking
**<u>*Masking functions for images and cubes*</u>**

**function mask_from_image, image, level, low** - create a mask from data *image*, setting to 1 if > *level*, to *low* (usually 0) if \<*level*

**function mask_from_image_nan, image, zero** - create a mask from data *image*, setting to 1 if  data value<>Nan. If *zero* = 0 or 1, set mask to Nan or 0 at Nan values.

**function mask_data, image, level, low** - masks data *image*, setting to *low* if < *level*

**function mask_data_median, image, level, low** - as above, but sets data image > *level* to median of *image*

**function mask_circle, data, x, y, r, v, rev** - masks *data* (image or cube) with circle center [*x,y*] radius *r*, set masked-out value to *v* (default 0). If *rev*<>0, reverse mask.

**function mask_set_nan_min, data, minvalue** - set *data* values to *minvalue* if value = Nan. If minvalue is zero, use the current minimum value. Equivalent to **cblank** function is *minvalue*=0

**function mask_cone, data, xc1, yc1, xc2, yc2, pa, beta, maskflag** - mask cone area over *data* (either image or cube), with equator [*xc1, yc1*], [*xc2, yc2*] (can be same coordinates for a point apex), centerline angle *pa* (from positive x-axis), internal full-angle *beta*. If *maskflag*=0, return the mask, if *maskflag*=1, return the masked input data.

**function mask_line, data, x1, y1, x2, y2, side** - creates a mask of *data* dimensions on one side of a line [*x1, y1*], [*x2, y2*]. *side* =0 for left, =1 for right side of line

### lib_velmap
**<u>*Velocity map (velmap) extension functions*</u>**

<a name="velmap_std_to_ext"></a>**function velmap_std_to_ext, velmapstd, r, cmin, vmethod, vzero, vx, vy** - convert standard QFitsView velmap *velmapstd* to extended form, *r*=instrumental resolution, *cmin*= minimum continuum value. Output is in extended velmap format - see below. Velocity zero is set by *vmethod* =

- 0 - median   
- 1 - average
- 2 - flux-weighted average
- 3 - manual (*vzero* value)
- 4 - pixel ([*vx,vy*] is set to zero)

**function velmap_vel_center, velmapstd, vmethod, vcenter, vx, vy** - returns the wavelength value from the standard *velmapstd* cube, using the methods as above

**function velmap_vel, velmapstd, vmethod, vcenter, vx, vy** - returns the velocity map from the standard *velmapstd* cube, using the methods as above

**function velmap_vel_set, velmapext, vmethod, vcenter, vx, vy** - Fix extended *velmap* velocity as per **velmap_std_to_ext** (re-do extended velmap cube)

**function velmap_rescale, velmapext, scale** - rescales extended *velmap* flux data (e.g. flux calib change)

<a name="velmap_fix"></a>**function velmap_fix, velmap, contlo, conthi, flo, fhi, vlo, vhi, wlo, whi, setvalue** - clean up *velmap* (either standard or extended form), setting values out of range to *setvalue*. Value ranges 

- *contlo, conthi* - continuum
- *flo, fhi* - flux
- *vlo, vhi* - wavelength
- *wlo, whi* - fwhm
- *setvalue* - value to set where spaxel is out of range (default 0/0)

**function velmap_extcorr, velmap, av, lambda** - extinction correct velocity map *velmap* at wavelength *lambda* (in nm), *av*=extinction A_V

**function velmap_extcorr_map, velmap, extmap, lambda** - as above, but *extmap* is a map of extinction values

<a name="velmap_fix_interp"></a>**function velmap_fix_interp, velmap, npix** - interpolate velmap *velmap* missing values, indicated by Nan in continuum layer (i.e. layer 1) (usually after **[velmap_fix](#velmap_fix)**). *npix* is interpolation width maximum

**function velmap_clean_map_wvt, velmap, map, nregion** - Clean up velmap *velmap* based on WVT *map* region number, setting region *nregion* pixels to NaN

**function velmap_mask, velmap** - set *velmap* to Nan where continuum=0

**procedure velmap_comps, velmapext, prefix, hmax** - Output *velmapext* components from *velmap*, to the current working directory. *prefix* (string) sets file names, terminated with 

- \_Flux - flux (layer 12)

- \_Flux\_Norm - normalized flux (range [0..1])
- \_Vel - velocity (layer 10)
- \_EW - equivalent width (layer 13)
- \_Sig - dispersion (layer 11)
- \_VelHist - velocity histogram. If *hmax* > 0, a 50 bin histogram of the velocity with range *hmax* -> *hmax*
- \_SigHist - - dispersion histogram. If *hmax* > 0, a 50 bin histogram of the velocity with range 0->*hmax*

**function velmap_from_profit, profit_data** - convert *profit_data* PROFIT cube format (see Riffel, R. A. 2010, Astrophys Space Sci, 327, 239, http://arxiv.org/abs/1002.1585) to standard velmap format.

**function velmap_derotate, cube, velmodel, lambdac** - Subtract velocity model *velmodel* from a data *cube*, where the velocity is determined from central wavelength *lambdac*. This shifts each spaxel spectrum by a wavelength amount calculated by the central wavelength and velocity model. This might be used if you have, say, a stellar rotation model and you want to apply it to gas emission lines.

### lib_chmap 
**<u>*Channel map functions*</u>**

<a name="chmap_create"></a>**function chmap_create, cube, lambda_cent, lambda_width, cutoff, width_factor, smooth** - make a channel map from the cube

- *lambda_cent* - estimate of central wavelength
- *lambda_width* - estimate of FWHM
- *threshold* - % of maximum for cutoff - default 0
- *width_factor* - wavelength widow (multiple of lambda_width) - default 2.5
- *smooth* - integer=boxcar, non-integer=gauss, 0=no smoothing - default 0

Returns a cube of channel maps, with axis 3 in velocity difference (km/s) from median. Spaxel values are FLUX (not flux density) in that channel

<a name="chmap_rebin"></a>**function chmap_rebin, cube, lnew, velwidth, sm, minval**- rebin channel maps in *cube* into *lnew* bins between velocities *v1* and *v2* (usually symmetric about 0, but not necessarily), with *sm* smoothing value, integer=boxcar, non-integer=gauss, 0=no smoothing, set output to NaN where < *minval*

<a name="chmap_comps"></a>**procedure chmap_comps, cube, dirout, fnameout** - splits channel map *cube* into components and writes images to folder *dirout*, named *fnameout* plus velocity (e.g. if *fnameout* ="pa_beta-450”, then output file name will be e.g.  "pa_beta-100.fits" etc.

**function chmap_from_velmap, velmapstd, cube_template, width, res** - create a channel map from a standard velmap *velmapstd*. The velmap is evaluated using *cube_template*, then the channel map is generated using the median velocity from the velmap with width multiplier (how far to extend the channels over the median FWHM) *width* (default 2.5). The FWHM is corrected by spectral resolution *res* (default 0).

### lib_pv
**<u>*Position Velocity Diagram functions*</u>**

**function pv_array, cube, ystart, wslit, nslit, lcent, lwidth** - create pv diagram from *cube* parallel to x axis, *ystart* - y pixel to start, *wslit* - slit width, *nslit* - number of slits, extract over range *lcent*-*lwidth* to *lcenter+lwidth*

**function pv_single, cube,  xc, yc, angle, width, lcent, vwidth, npix, contflag** - extract single PV plot at *xc*/*yc*/*angle*/*width* - centerered on *lcent*. *vwidth* - velocity width around *lcent*, rebinned in velocity to *npix* channels. *contflag* =1 subtract continuum (flux) =2 divide continuum (effectively the same as equivalent width) =0 don’t remove continuum

**function pv_ratio, cube,  xc, yc, angle, width, lcent1, lcent2, vwidth, npix** - create PV diagram as above for ratio of 2 lines *lcent1*, *lcent2*

**function pv_meddev, image** - divide *image* by median along x axis (useful for EW for PV diagrams)

### lib_wvt
**<u>*Weighted Voronoi Tesselation functions*</u>**

**function wvt_cube, cube, sn_target** - make WVT *cube* using noise in each spaxel, to *sn_target*. Bad pixels where S/N is > 10x brightest pixel S/N

<a name="wvt_cube_mask"></a>**function wvt_cube_mask, cube, l1, l2, mask, cutoff, sn1, sn2** -  make WVT cube using 2 S/N ratios, inside and outside *mask*. Returns WVT applied to *cube*.

- *l1, l2* - wavelength range to use for signal and noise determination (“quiet” part of spectrum with no emission lines)
- *mask* - if 2D mask, use this. If *mask*=0, use *cutoff* to determine mask
- *cutoff* - percentage of peak maximum for mask level
- *sn1*, *sn2* - S/N ratios for inside/outside mask. If *sn2*=0, just use *sn1* over whole cube

**function wvt_sn_mask, cube, l1, l2, mask, cutoff, sn1, sn2** - as for **[wvt_cube_mask](#wvt_cube_mask)**, except returns WVT image data - layers:

1. Signal
2. Noise
3.  S/N
4.  Mask
5.  Signal binned
6.  Signal bin map
7.  Bin density (1=maximum - smallest bins , 0=minimum - biggest bins)

**function wvt_build_from_map_cube, cube, wvtmap, prntflag** - make WVT cube from *cube* and *wvtmap*. If *prntflag* =1 print diagnostic every 100 regions

**function wvt_build_from_map_image, image, wvtmap** - make WVT image from *image* and *wvtmap*.

**function wvt_velmap, velmap, layer, sn** - make WVT velmap from standard or extended velmap *velmap*, *layer* is either 0=continuum, 1=flux, *sn*=S/N target

**function wvt_density, wvtmap** - make map of region density, i.e. 1/# of pixels in region. *wvtmap* is WVT with /map flag.

**function wvt_cube_to_specarray, cube, wvtmap, normflag, prntflag** - convert *cube* inbuff to spectrum array, using *wvtmap* regions. If *normflag* = 1, divide each spectrum in array by the first one. If *prntflag* = 1, print running diagnostics

**function wvt_specarray_to_cube, image, wvtmap** - reverse of **wvt_cube_to_specarray**

### lib_general

***<u>Miscellaneous functions</u>***

**function indexreform, index, xsize, ysize, zsize** - returns 3D co-ords from 1D *index*, given dimensions *xsize, ysize, zsize*. Values returned as array.

**function lognan, data** - set log of *data*, setting zero and Nan values to Nan

**function clipnan, data, low, high** - set values outside range [*low..high*] to Nan

**function axiscentroids, image, axis** - returns centroids of each image row/column, row-*axis*=1, column-*axis*=2; used e.g. for finding centroids of pv diagram

**function histogram_bin, data, low, high, bin, normflag** - create a histogram from inbuff data (any dimensions), from _low_ to _high_ values in *bin* bins. Histogram is normalised if *normflag* =1. Output x-axis values are set to range.

**function profile_export, data, scale1, scale2, scale3, offset** - Export 1D profiles from *data* with up to 3 separate scales, e.g. arcsec, pc, Re plus the pixel scale, _offset_=1 offsets by 1/2 a pixel (e.g. for log scale plot) (default 0). Scales default to 1 if not given.

**function butterworth_filter, order, cutoff, size** - Create a Butterworth filter for order _order_ for a square of sides _size_, with _cutoff_ Nyquist frequency.

**function interp, data, x, x1, x2**  - linearly interpret over *data* at position *x* over *x1*-*x2*.
Used to image_interp_flags and cube_interp_flags

### lib_astronomy

**<u>*General astronomy functions, implemented from IDL.*</u>**

**function G** - gravitational constant (MKS)

**function Msun** - Mass of the sun in kg

**function Pc** - 1 Parsec in meters

**function airtovac, wave** - Convert air wavelengths to vacuum wavelengths, *wave* in Å

**function planck, wave, temp** - Calculates the Planck function in units of ergs/cm2/s/Å. *wave* in Å, *temp* in degrees K.

**function coordstring, ra, dec, rad** - create a nice string of celestial coordinates. *ra* and *dec* should be given in radians; if *rad*>0, convert to degrees. Example:
`print coordstring(9.1234,5.678,1)`
`RA = 0h 36m 29.62s, DEC = + 5d 40' 40.8"`

### lib_astro_general

**<u>*General astrophysics functions.*</u>** 

All the "lib\_astro\_*.dpuser" functions are executed from the "lib\_astro.dpuser" script.

**function redshift_data, data, z, rdflag, nanflag, smth, fcons** - redshift data by *z*, assuming last axis is wavelength. WCS values set. If *rdflag* = 1, redisperse shifted data to same wavelength range as input. If *nanflag*=1, set nans in output at same pixels as in input. If *smth*>0, smooth by no. of pixels. If *fcons*=1, conserve total values.

**function bb_make,t, l1, l2, npix, wlflag** - make black-body function at temparture *t*, wavelength range *l1* to *l2*, number of pixels *n*, *wlflag*=0, wavelength in A, =1=> nm, =2 => um

**function bb_make_log, t, l1, l2, npix, scale, cutof**f - make bb at temp *t* over log wavelength [*l1,l2*] (in log meters), creating spectrum length *npix*, multiply wavelengths by *scale* (to convert to e.g. nm), set result to Nan where below *cutoff*. 

**function bb_div, spectrum, temp** - divide *spectrum* by black-body at temperature *t*

**function extinction_calc, f1, f2, l1, l2, rat, galext, s, flmin1, flmin2**- create an extinction map from 2 emission line maps - *f1*, *f2*  are flux maps, *l1*, *l2*=wavelengths, *rat*=expected flux ratio, *galext*=galactic extinction, *smth*= smooth pixels, *flmin1*, *flmin2*=minimum flux value for each map - calculates the extinction constant (CCM laws for IR and optical)

**function extinction_correct, cube, av** - correct *cube* for extinction (*av* =single value for extinction) - wavelength from axis 3

**function extinction_correct_map, cube, av** - correct *cube* for extinction (*av* =extinction map) - wavelength from axis 3.

**function extinction_correct_lambda, data, av, lambda** - correct value/image for extinction *av* at wavelength *lambda*: can be used on value or image

### lib_astro_mapping
**<u>*Astronomy functions (mapping and excitation diagrams)*</u>**

**function map_compare_diagram, image1, image2, min1, max1, min2, max2, nbin, lgaxesflag** - Map diagram density plot. *image1*, *image2* - value maps, x and y axes. *min1/2*, *max1/2* - min and maximum values for axes 1/2. *nbin* - no of bins on each axis. *lgaxesflag* - 1=plot in log space (min,max must be in log values). Generates e.g. BPT diagrams

**function map_compare_pos, image1, image2, image3, image4, x, y, boxsize** - get 2 sets of map ratios (*image1*/*image2*, *image3*/*image4*)at position [*x,y*], averaged over *boxsize* x *boxsize* pixels (e.g. excitation ratios at feature position)

**function map_basis_distance, basex0, basey0, basex100, basey100, x1, x2, y1, y2, size** - creates an image (dimensions *size*, limits (*x1, y1*), (*x2, y2*)) of distance from basis points 0 to 100% [*basex0, basey0*] to [*basex100, basey100*] - for use in AGN mixing ratios for contour values.

**function map_compare_basis, image1, image2, basex0, basey0, basex100, basey100, lgaxesflag** - plots basis distance (AGN mixing ratio) from basis points [*basex0, basey0*] to [*basex100, basey100*] . *lgaxesflag* - 1=take log of *image1*, *image2* before calculation

**function map_regime_ir, image1, image2, a1, a2, a3, b1, b2** - create position excitation map. If a1=0, use the standard infrared Riffel 2013 excitation regimes. *image1* is H_2/Br_gamma, *image2* is [Fe II]/Pa_beta. Both in log values. Output values at each spaxel are SF=1, AGN=2, LINER=3, TO1=4, TO2=5

**function map_regime_optical, image1, image2, typeflag**- create position excitation map for optical line ratios (*image1* and *image2*) from Kewley et al. 2006 regimes. *typeflag* = 1 ([N II]/H_alpha diagram), =2 ([S II]/H_alpha diagram), =3 ([O I]/H_alpha diagram). Returns 1=SF, 2=Seyfert, 3=LINER, 4=Composite

### lib_astro_spectrum
**<u>*Astronomy functions (spectrum)*</u>**

**function spec_fluxdens, spectrum, l1, l2, prflag** - flux density (counts/nm) for *spectrum* between *l1* and *l2* <u>wavelength</u>; returns a single value. If *prflag*<>0, print results as well.

**function spec_sn, spectrum, l1, l2** - Compute S/N for *spectrum* over wavelength region *l1*-*l2*, using median rather than average, as more robust.

**function spec_wave_to_vel, spectrum, lambda -** Convert *spectrum* wavelength axis to velocity, with zero velocity wavelength *lambda* .

### lib_astro_image
**<u>*Astronomy functions (image)*</u>** 

**function img_aphot_annular, image, xcen, ycen, r, ib, ob** - aperture photometry on *image*,centered on [*xcen,ycen*]; aperture *r*, background annulus from *ib* to *ob* (inner to outer boundary). If *ib* and *ob* are zero, set to r and 2*r

**function img_apphot_simple, image, xcen, ycen, r, pixsize, scale** - simple aperture photometry on image,centered on  [*xcen,ycen*]; aperture *r*.

**function img_flux_to_mag, image, zpm, ssize, zpflag** - convert flux *image*  to mag, *zpm* is either zero-point magnitude (*zpflag*=0) or zero-magnitude flux (*zpflag*=1 - default). *ssize* = pixel size in arcsec to convert to mag/arcsec^2 (default 1). This can also work for a single number, in which case *ssize* is set to 1.

**function img_convert_filters, image1, image2, coeffs12, coeffs21, tol** - Convert images in 2 filters to another filter set using the methodology of ﻿Holtzman, J. A., Burrows, C. J., Casertano, S., et al. 1995, PASP, 107, 1065. The coefficients (*coeff12* and *coeffs21*) for WFPC2/WFC3 are from Holtzmann, for ACS from﻿ Sirianni, et al. 2005, PASP, 117, 1049. The tolerance for the result is where the maximum magnitude change on each iteration is less than *tol*. 

Example - for conversion from ACS (WFC) F606W and F814W images to *V−I* colour image, the coefficient sets (Table 22 of Sirianni et al.) are:
*coeffs12* - [26.325, 0.236, 0.000]
*coeffs21* - [25.485, -0.002, 0.000]

The function returns a cube with 3 layers, (1) filter result 1 e.g. *V* (2) filter result 2 e.g *I* (3) difference e.g. *V-I*, as well as printing iteration and result diagnostics.

**function img_SFD_dust_pos, dustimgn, dustimgs, lat, long** - returns galactic dust extinction value, as per Schlegel et al. 1998, ApJ, 500, 525. *lat* and *long* are galactic co-ordinates. *dustimgn* and *dustimgs* are the dust mappings, north and south. The standard for *E(B-V)* is the "SFD_dust_4096" maps, but others can be used - these are all downloadable from https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/EWCNL5.

**function img_SF_dust_pos, dustimgn, dustimgs, lat, long** - as above but corrects the *E(B-V)* as per Schlafly, E. F., & Finkbeiner, D. P. 2011, Astrophys J, 737, 103.

### lib_astro_cube
**<u>*Astronomy functions (cube)*</u>** 

**function cube_apspec, cube, ox, oy, or, bx, by, br, br2, mask** - Get star spectrum from *cube* withusing circular aperture and background, plus a mask circle. 

- *ox, oy, or* - center and radius of aperture
- *bx, by, br* - center and radius of mask, if not required then bx = 0
- *br2* - inner radius of mask annulus (default 0). This annulus allows the background to be centered the same as the aperture (in which case *br2* > *or*)
- *mask* - any other mask required (2D fits) (default none)

Output is spectrum of aperture less average of background (with mask) - values < 0 are set to Nan.

**function cube_fluxdens, cube, l1, l2, prflag/function spec_fluxdens, spectrum, l1, l2, prflag** - flux density (counts/nm) between *l1* and *l2* wavelength; returns image (cube\_) or single (spec\_) value. If *prflag*<>0, print results as well

**function cube_sky_rem, cube, bckgnd_lvl** - removes skylines from *cube*. Takes background pixels as those with median value below *bkgnd_lvl*

**function cube_sl_clean, cube, skyline_list, width** - removes skylines from *cube* using skyline_linelist (array of wavelengths) - interpolated over wavelength ± *width*

<a name="cube_clean_bp_fix"></a>**function cube_clean_bp_fix, cube, bp_cube** - cleans *cube* based on bad pixel cube *bp_cube* using **dpixapply** over x image slices

**function cube_clean_bp, cube, threshold** - create bad pixel cube using *threshold* scanning over wavelength slices.

**function cube_clean_bp_limits, cube, ll, ul** - create bad pixel cube from *cube* for input to **[clean_cube_bp_fix](#clean_cube_bp_fix)**, flagging pixels below *ll* and above *ul* values

### lib_astro_kinematics

**function velmodel_plummer, x0, y0, M0, Re, psi0, inc, xsize, ysize, pixscale, angscale** - create a Plummer kinematic rotation model velocity field in km/s. This can be applied to a cube by e.g. **[cube_velocity_correct](#cube_velocity_correct)**.
*x0, y0* - centre pixel position
*M0* - enclosed mass (units of Msun)
*Re* - length scale (pc)
*psi0* - line of nodes (default 0 - counterclockwise from +ve X axis)
*i* - inclination of disk (0-90, 0 = face-on - default)
*xsize, ysize* - image size (by default, twice x0, y0)
*pixscale* - pixel scale in arcsec (default 1)
*angscale* - angular scale in pc/arcsec (default 1)

**function velmodel_hernquist, x0, y0, M0, Re, psi0, inc, xsize, ysize, pixscale, angscale** - as for **velmodel_plummer** for Hernquist kinematic rotation model.

**function velmodel_geom, x0, y0, psi0, inc, xsize, ysize** - computes the geometric components of rotational models, with parameters as above. Returns data [*xsize, ysize*, 2] - first layer is rotational component, second layer is radial component. These are independant of the particular model form.

# Standard Procedures

## Velocity maps
* Create QFitsView **velmap** with wavelength, fwhm estimate
* Examine velmap for continuum, height, wavelength and fwhm “sensible” ranges
* Use **[velmap_fix](#velmap_fix)** to clean up velmap, entering "good" ranges for each fit component. This sets out-of-range spaxels to Nan.
* Use **[velmap_fix_interp](#velmap_fix_interp)** to interpolate over NaN values (if required)
* Use **[velmap_std_to_ext](#velmap_std_to_ext)** to create extended velmap format

**<u>*Standard VELMAP Format*</u>**

1. Continuum
2. Peak height above continuum
3. Wavelength
4. FWHM
5. e_Continuum
6. e_Peak
7. e_Wavelength
8. e_FWHM
9. Chi-squared

**<u>*Extended VELMAP Format*</u>**

1. Continuum
2. Peak height above continuum
3. Wavelength
4. FWHM
5. e_Continuum
6. e_Peak
7. e_Wavelength
8. e_FWHM
9. Chi-squared
10. Velocity (zero-point calculated/set by *method* parameter in the **[velmap_std_to_ext](#velmap_std_to_ext)** function)
11. Dispersion (sigma) velocity, corrected for spectral resolution
12. Flux (Peak\*FWHM*1.0699)
13. Equivalent width (flux/continuum)
14. Total support (order + turbulence) (√)
15. Order vs turbulence (|V/σ|)

## Channel maps

- Create basic channel map using **[chmap_create](#chmap_create)** (usually do not smooth)
- Rebin to required # of channels (e.g. 9 or 16) using **[chmap_rebin](#chmap_rebin)** (smoothing if required)
- Output individual channel maps using **[chmap_comps](#chmap_comps)**
