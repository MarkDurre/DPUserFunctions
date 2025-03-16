# QFitsView DPUser
## Setting Up DPUser Functions Library
1. Place all libraries (files with “lib\_\*.DPUser” name) in a convenient location  e.g. “*my\_code\_location*/DPUserlib/Functions" (subsituting for *my\_code\_location* e.g. "/Users/*username*/Programs").

2. Place “startup.DPUser” in e.g. “*my\_code\_location*/DPUserlib". This assumes that all the library functions are located in the "Functions" sub-folder and runs the "lib\_all.DPUser" script. This file must be modified for your own requirements; it also sets the *DPUser\_DIR* environment variable that can be accessed by other scripts, using the **getenv** function in QFitsView.

    ```
    //Example startup.DPUser
    //This line must be modified (change my\_code\_location) for individual users
    setenv "DPUser_DIR", "my_code_location/DPUserlib"
    //DPUserdir is used by libraries to call other libraries if required
    DPUserdir=getenv("DPUser_DIR")
    print "Running General Functions : DPUser Directory - "+DPUserdir
    run DPUserdir+"/Functions/lib_all.DPUser"
    print "Finished General Functions - "+DPUserdir
    //You can put your own startup DPUser code here
    ```

3. As of version 4.1 of QFitsView, the location of the "startup.DPUser" file can be specified in the QFitsView menu "QFitsView > Preferences > Paths > DPUser Path" (Mac) or "Options > Preferences > Paths > DPUser Path" (Windows and Linux). 

    For version previous to 4.1, create a symbolic link in the root directory to this folder. For macOS 10.15+ use the `/etc/synthetic.conf` symbolic links method - reference [here](https://stackoverflow.com/questions/58396821/what-is-the-proper-way-to-create-a-root-sym-link-in-catalina)). This link must be called "DPUserlib".

4. When QFitsView starts, it automatically looks for and executes the script file "startup.DPUser" in the DPUser Path specified above (for versions before 4.1, it looks for "/DPUserlib/startup.DPUser") i.e. it runs the script set up above. This runs all libraries to make the functions available to QFitsView - you will see a whole bunch of “Stored function…” and “Stored procedure…” plus “Finished General Functions” text lines on the DPUser area in QFitsView.

5. If the above folder conventions are not used, the following files must be modified:
    - *my\_program\_location*/DPUserlib/startup.DPUser
    - *my\_program\_location*/DPUserlib/Functions/lib\_all.DPUser

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

## Editing DPUser Code with BBEdit

The main documentation for DPUser is through [this link](https://www.mpe.mpg.de/~ott/DPUser/). DPUser code can be edited with QFitsView *DPUser > Script Editor*. It can also be edited by various external text editors; e.g. for Macs, **BBEdit** (as well as using the internal script editor). To facilitate this, a language module has been implemented for BBEdit - "DPUser.plist", with the following highlighting features: 

- Syntax - both structural commands - e.g. "if", "else" etc. and internal DPUser functions/procedures. As new functions/procedures are implemented in QFitsView, the "BBLMPredefinedNameList" array must be updated.
- Comments (both for "/\*..\*/" and "//").
- Strings ("..")
- Function and procedure prefixes

The file is installed in BBEdit's language module directory - by default "/Users/*username*/Library/Application Support/BBEdit/Language Modules".

Note that after editing code in an external text editor, the script must be executed again with QFitsView to make functions/procedures available.

##  Parameter data types
In the following description of functions/procedures, 

* FITS or buffer input - *cube* is 3D, *image* is 2D, *spectrum* is 1D, *data* is 1, 2 or 3D, *mef* is multi-extension FITS data (created using the **list** function for up to 10 data buffers).
* Pixel co-ordinates
   * p, p1, p2,… (general pixels co-ordinates)
   * x, x1, x2,…, y, y1, y2,…, z, z1, z2,… (for x, y or z axes)
   * pixel/wavelength masking - an array of 2n values [x1,x2,x3,x4...] - paired pixel numbers or wavelengths (pairs p1..p2,p3..p4 etc.) for spectral masking
* WCS co-ordinates
   * wcs - for axis set [CRPIXn, CRVALn, CDELTn]
  * w, w1, w2,…. (for individual coordinates)
  * l, l1, l2, …. (wavelengths)
* Maps
   * velmap - velocity map cube (either standard *velmapstd* or extended *velmapext* format)
  * wvtmap - weighted Voronoi tessellation map (region numbers)
* NaN - "not a number", defined in DPUser as 0/0.

## Libraries
### lib\_all

Runs all the following libraries; this just consists of script lines to execute other scripts in the "Functions" sub-folder. Other scripts can be executed by adding appropriate lines, e.g.

`@SomeFolder/SomeScript.DPUser`

If full path to script is not given, it is assumed to be relative to the DPUserlib folder defined above. The DPUser command *@* means execute this line as a script.

### lib\_header

***<u>Enhanced functions on FITS headers</u>***

**function hdr\_get\_fits\_key, data, key** - replaces **getfitskey** function with a check that *key* exists. If it is not found, then returns blank. This is required since QFitsView function _getfitskey_ returns an error if the key is not found.

**function hdr\_get\_lines, data** - returns string array of all lines in *data* header (rather than single string that function "header" does)

**function hdr\_check\_prefix, data, key** - checks *data* header for presence of *key* at start of line. If found, returns the line number - this only matches the first found. If not found, this returns 0. This function can be used to check if the *key* exists before using _getfitskey_ and potentially getting an error.

**function hdr\_get\_all\_fits\_keys, data, prefix** - returns a string array of all keys in *data* header that match prefix, e.g. if *prefix* = "ZPT\_" then the array will contain ["ZPT\_001", "ZPT\_002"...]. If *prefix* is not given, it is set to blank, i.e. all keys are matched.

**function hdr\_get\_all\_fits\_key\_values, data, prefix, type** - as above, but returns the values of those keys. If *type* = 0 (default), returns a string array. If *type* = 1, returns a numerical array; string keys are returned as NaN, i.e. 0/0. 

**function hdr\_read\_from\_text, hdrdata, prefix** - read header keys from *hdrdata* (stringarray) and add to FITS file. Default *prefix* for header keys is "#". Normally the header fits file will use DPUser _copyheader_ to transfer to data fits file, but this is for some web file downloads.

### lib\_wcs

**<u>*Transform to and from World Coordinate Systems and Pixels*</u>**

**function get\_WCS\_data, data, axis** - return array [CRVAL, CRPIX, CDELT] for *axis* (1,2 or 3) of *data* 

**procedure set\_WCS\_data, data, wcs, axis** - sets WCS values for *data* for *axis* (1,2 or 3)

**function cvt\_pixel\_WCS, pix, crpix, crval, cdelt** - convert pixel number *pix* to WCS coordinates using *crpix, crval, cdelt*, asssuming linear projection (cartesian). The new DPUser _worldpos_ function does this including the correct projections.

**function cvt\_WCS\_pixel, value, crpix, crval, cdelt** - convert WCS coordinate *value* to pixel using *crpix, crval, cdelt*, asssuming linear projection (cartesian).  The DPUser _pixpos_ function does this including the correct projections.

**function cvt\_WCS\_pixel\_data, data, value, axis** - converts pixel *value* to WCS for *data* for *axis*

**function cvt\_pixel\_WCS\_data, data, value, axis** - reverse of *cvt\_WCS\_pixel\_data*

**function WCS\_range, data, p1, p2, axis, prntflag** - returns WCS coordinates as range [*w1,w2*] from pixel values [*p1,p2*] for axis on *data*, if *prntflag*=1, then print range

**function pixel\_range, data, w1, w2, axis, prntflag**  - returns pixel values as range [*p1,p2*] from WCS co-ordinates [*w1,w2*] for axis on *data*, if *prntflag*=1, then print range

**function set\_WCS\_default, data** - checks *data* has minimal WCS keys set (to 1 by default).

<a name="get\_WCS\_image"></a>**function get\_WCS\_image, image** - Get axis 1 and 2 WCS data for *image* (can be cube) and calculates CDELT and rotation angle from CD keys; result is array of key values [CRPIX1, CRVAL1, CD1\_1, CD1\_2, CRPIX2, CRVAL2, CD2\_1, CD2\_2, CDELT1, CDELT2, CROTA2]

**function set\_WCS\_image\_scale, image, wcs2d, xscale, yscale** - Rescales (e.g. for non-integer re-binning) using *xscale* and *yscale* and sets WCS data for *image* (or cube), including CD keys - removes CDELT and CROTA2 keywords. Input *wcs2d* is format as for **[get\_WCS\_image](#get\_WCS\_image)**.

**function get\_WCS\_values, data** - create WCS array [1,cv,cd] from dispersion *data*, calculated from first/last values and number of elements.

**function WCS\_cdelt\_cd, cdelt1, cdelt2, rotang** - converts WCS x,y pixel sizes (*cdelt1, cdelt2*) and rotation angle (*rotang*) to CD matrix values. Returns vector [CD11, CD12, CD21, CD22].

**procedure set\_cd\_keys, data, cdkeys** - sets CD keys in FITS header of *data* from vector *cdkeys* (in same format as **WCS\_cdelt\_cd** function) and deletes the CDELT1/2 keys.

**function WCS\_shift\_pix, data, xshift, yshift, sec** - shift image or cube astrometry by altering the CRPIX1,2 values (useful to align images). *xshift, yshift* - amount to shift in axis 1 and 2 respectively, *sec* (default 0), if =0, shifts in pixels, else in seconds of arc. Note for seconds of arc shift, you must mutiply RA seconds of time by 15.

**procedure copy\_WCS\_data, data1, data2, axis1, axis2** - copy WCS values from *data1*/*axis1* to *data2*/*axis2*.

**function WCS_PC_To_CD, inbuff** - converts PC/CDELT WCS parameters  to CD matrix (preferred).


### lib\_cube
**<u>*Data cube functions*</u>**

**function cube\_trim\_xy, cube, x1, x2, y1, y2** - sets cube to zero for x<*x1*, x>*x2*, y<*y1*, y>*y2* (_cblank_ cube first)

**function cube\_trim\_wl, cube, w1, w2, value, trimflag** - sets *cube* to *value* (usually zero) for l < *l1*, l > *l2* in axis 3 using WCS (_cblank_ cube first). If *trimflag*=1, then truncate the cube outside the wavelength range.

**function cube\_spectrum\_mask, cube, mask, level** - mask *cube* on spectral wavelength with *mask* (pixel pairs), set masked pixels to *level*

**function cube\_clip, cube, lvl, thresh, mask** - clips *cube* <0 and > *lvl* in image (x-y) plane, then does _dpixapply_ using threshold, cube is spectrally masked

**function cube\_clip\_y, cube, lvl, thresh** - as above, but in the x-z plane    

**function cube\_interp\_z, cube, x1, x2, y1, y2, z1, z2** - interpolate *cube* in image plane over rectangle [*x1:x2,y1:y2*] in each of wavelength plane [*z1:z2*]

<a name="cube\_interp\_x"></a>**function cube\_interp\_x, cube, x1, x2, y1, y2, z1,z2** - interpolate *cube* in image plane over rectangle [*y1:y2,z1:z2*] in each of spatial range [*x1:x2*]

<a name="cube\_interp\_y"></a>**function cube\_interp\_y, cube, x1, x2, y1, y2, z1,z2** - interpolate *cube* in image plane over rectangle [*x1:x2,z1:z2*] in each of spatial range [*y1:y2*]

<a name="cube\_interp\_xy"></a>**function cube\_interp\_xy, cube, x1, x2, y1, y2, z1,z2** - as above, but interpolate over wavelength [*z1:z2*] in the xy plane

**function cube\_set\_value, cube, x1, x2, y1, y2, z, xv, yv** - set rectangle [*x1:x2,y1:y2*] at image plane *z* to value at [*xv, yv*]

**function cube\_pixfix\_xy, cube, pixfixdata, n** - fix cube using **cube\_interp\_xy** function, n sets, *pixfixdata* is  *n* x 6 array [x1,x2,y1,y2,z1,z2].

**function cube\_single\_pixel\_fix, cube, x, y** - do **cube\_interp\_xy** for all z axis for single spaxel

**procedure cube\_bit\_nan, cube,x,y**  - set spaxel [*x*,*y*] to 0/0 along whole *cube* z axis

**function cube\_clean\_dpix, cube, scale** - Clean *cube*  by dpixcreate/apply, the threshold for dpixcreate is set from maximum of median image divided by *scale*.

**function cube\_resize\_center, cube, xcent, ycent, xsize, ysize, subpix** - resize *cube*  to size *xsize*,*ysize* and center on pixel [*xcen,ycen*]. By default, *ycent=xcent*, *xsize, ysize*= *xcent\**2. If *subpix*=1, perform sub-pixel shifting on the cube, else shift by integer pixels (default)

**function cube\_shift\_xy, cube, xshift, yshift** - sub-pixel shifts *cube* by [*xshift, yshift*]

**function cube\_shift\_z, cube, zshift** - sub-pixel shifts *cube* by [*zshift*] in the z (wavelength) axis.

**function cube\_redisp, cube, disp\_old, disp\_new, prnt** - redisperses *cube* (axis 3) to new from *disp\_old* to *disp\_new* dispersion spectra by interpolation. If *prnt* <> 0, print diagnostics.

**function cube\_redisp, cube, spectrum** - Simple redispersion of *cube* spectral axis to *spectrum* dispersion.

**function cube\_symm\_flip, cube, lambda, width, part** - symmetrically flip *cube* about wavelength *lambda*, *part* =0 (left) or 1 (right), trims cube to *lambda*+-*width*

**function cube\_rotate, cube, xcen, ycen, rot\_angle, pixscale** - rotate *cube* on center [*xcen*,*ycen*] by *rot\_angle*, setting *pixscale* in arcsec/pixel. This is used because the **rotate** function does not work on cubes.

**function cube\_centroids, cube** - get centroids at each wavelength layer (z axis). Returns a FITS array of dimensions naxis3(*cube*) x 2, with x and y centroids at each pixel layer. The wavelength WCS is set.

<a name="cube\_centroids\_gauss"></a>**function cube\_centroids\_gauss, cube, xe, ye, we, mask** - get centroid at each pixel layer, with estimated center at [*xe,ye*] over fitting window *we*. The spectrum is masked by *mask* (if not zero or not entered).

**function cube\_centroid\_gauss\_align, cube, xc, yc, xe, ye, we, mask** - align *cube* centroids at each pixel layer, using [*xe, ye, we, mask*] are the estimate parameters of the peak (as for **[cube\_centroids\_gauss](#cube\_centroids\_gauss)**), with the centroids aligned to [*xc, yc*].

**function cube\_cont\_slope, cube, mask** - returns image with continuum slope at each spaxel of *cube*, masked by wavelength pairs *mask*

**function cube\_spectrum\_add, cube, spectrum, x1, x2, y1, y2** - adds *spectrum* to *cube* for each spaxel. If any of *x1, x2, y1, y2* are set to zero, then perform the action over the pixel range. By default (not given), these are set to zero. To subtract *spectrum*, add by negative.

**function cube\_spectrum\_multiply, cube, spectrum, x1, x2, y1, y2** - multiply *cube* by *spectrum* as for **cube\_spectrum\_add**. To divide by *spectrum*, multiply by inverse.

**function cube\_set\_pixlayers, cube, pixl, p1, p2** - set *cube* layers [*p1, p2*] to the values for layer *pixl*.

**function cube\_wavelength\_correct, cube, correction** - corrects the wavelength solution at each spaxel by shifting the spectrum.  *correction* is an image of the same dimensions as the *cube* x and y axes and is in same units as the spectral axis of *cube*. 

<a name="cube\_velocity\_correct"></a>**function cube\_velocity\_correct, cube, velmodel** - applies a velocity field model *velmodel* (in km/s) to each spaxel of *cube*. Each spaxel is re-dispersed by intepolation.

**function cube\_to\_2d, cube** - Convert data *cube* to 2d apertures for IRAF. Returns a 2D array with spectrum on the x-axis and all spaxels on the y axis.

**function cube\_set\_flags\_NaN, cube, layer** - set up flags image for cube\_interp\_flags, from a data *cube* (e.g. a velmap) from *layer*. This retuens an image with same dimensions as x and y axes as the cube, with 1 where pixel in NaN, 0 else.

**function cube\_interp\_flags, cube, flags, xi1, xi2, yi1, yi2, dmax** - interpolate over pixels in *cube* where *flags* is set to 1, 0 = good values to use for interpolation. [*xi1:xi2, yi1:yi2*] is region to interpolate (*xi1* = 0 - do whole area). *dmax* is maximum distance from “good” pixels (by default = 1). *flags* can be generated from **cube\_set\_flags\_nan**.

**function cube\_deslope, cube, mask, wlflag** - deslope *cube* for each spectrum using **spectrum\_deslope**. *wlflag* = 1 if mask values in wavelength

**function cube\_clean\_pixels, cube, layer, npix** - Remove singleton pixels surrounded by NaN’s, opposite of **cube\_interp\_flags**, used to clean up boundaries etc. *npix* is max number of good pixels around each pixel before blanking.

**function cube\_radial\_spectrum, cube, xc, yc, rstep, nstep, ann** - Radial spectra of *cube*, centered [*xc,yc*] radial steps *rstep*, number of steps *nstep*. If *ann=1*, output annular spectra

**function cube\_from\_image\_spectrum, image, spectrum** - Creates a cube from an *image* and *spectrum*. Wavelength axis of cube is spectrum scaled by image value.

**function cube\_rebinxy, cube, xscale, yscale, kernel** - Rebin *cube* or image pixel scaling in x and y directions by *xscale*, *yscale*. Uses the _interpolate_ DPUser function with kernel *kernel*. Note this function DOES NOT handle the WCS co-ordinates scaling; use **get\_WCS\_cube** and **set\_WCS\_cube\_scale** functions.

**function cube\_rebinfrac, cube, xscale, yscale** - Rebins *cube* (image) to *xscale*, *yscale* using fractional binning. Note comments about WCS values as above.

**function cube\_rebinx, cube, xscale**  - Rebin *cube* (image) pixel scaling in x direction ONLY. Uses the **interpol** DPUser function (quicker than **interpolate**). Note comments about WCS values as above.

**function cube\_combine\_avg, mef, omit** - combine cubes by averaging. *mef* is a multi-extension fits list, *omit* is a value to reject in the averaging. If *omit* is not provided or set to 0/0, then no rejection is done. *mef* is created from mutiple cubes (up to 10 - this is a DPUser limit) by e.g. 
`mef = list(buffer1, buffer2, buffer3, ....)`

**function cube\_combine\_median, mef, omit** - median combine cubes as for **cube\_combine\_avg**.

<a name="cube\_apply\_snr"></a>**function cube\_apply\_snr, cube, snr, snrscale** - adds noise to a signal *cube* using signal-to-noise ratio *snr*, multiplied by factor *snrscale* (default 1). *snr* can be a single value, a spectrum which is applied at each spaxel, an image where a single value is applied at each spaxel or a cube with different values at each pixel. The noise is a random gaussian value with standard deviation of the applicable SNR.

**function cube\_spectrum\_sub\_scale, cube, spectrum, l1 ,l2**  - Subtract scaled *spectrum* from *cube*, for spectral model continuum removal, *l1, l2* is wavelength range for scaling. The *spectrum* must have the same wavelength axis size as *cube*

<a name="cube\_flux\_ul"></a>**function cube\_flux\_ul, inbuff, lambda, width, deslope, ignore, flag** -  Calculate flux upper limits. Returns SD and (optionally) average for *cube* over wavelength range defined by central wavelength *lambda* and window *width*. This function assumes spectrum is in last axis and has correct WCS. *deslope*  - order of polynomial fit to slope, if zero, do not deslope (default). *ignore* - value to ignore (e.g. 0 or NaN) - default NaN, i.e. don't ignore any number. *flag* - if 0 (default) return only SD, if 1, return [SD, AVG]

**function cube\_centroids, cube, pos, radius** - Find centroid along z axis for *cube*, initial estimate at *pos* [x,y] co-ordinates. Allow search radius about pos. Returns [x,y] arrow of positions at each wavelength. Uses the DPUser function _centroids_.

**function cube\_3d\_to\_2d, cube, axis** - Convert *cube* from 3d to 2d image. Scanned along either x (*axis*=1 - default) or y (*axis*=2)

### lib\_image
**<u>*Image functions*</u>**

**function image\_erodenan, image** - erode *image*, pixels set to NaN if any neighbour is NaN.

**function image\_smooth, image, smooth** - smooth *image* that has NaN values - *smooth* integer=boxcar, non-integer=gaussian. DPUser functions would fail with NaN values

**function image\_interp\_x, image, x1, x2, y1, y2** - as for **[cube\_interp\_xy](#cube\_interp\_xy)**, but for single *image* .

**function image\_interp\_y, image, x1, x2, y1, y2** - as for **[cube\_interp\_x](#cube\_interp\_x)**, but for single *image* .

**function image\_interp\_xy, image, x1, x2, y1, y2** - as for **[cube\_interp\_y](#cube\_interp\_y)**, but for single *image* .

**function image\_from\_profile, profile, xp, yp, xc, yc** - create 2D image from 1D *profile*, size of output image is *xp* x *yp* , [*xc, yc*] - center of rebuilt profile

**function image\_bfilter, image, order, cutoff** - Butterworth filter an *image*, assume square image, filter order *=order*, *cutoff* =Nyquist cutoff (0-1)

**function image\_enclosed\_flux, image, xc, yc, r, smth** - Get enclosed flux within radius *r* from [*xc, yc*] (pixels). If *smth*>0, Gaussian smooth the output 

**function image\_avg, image, x, y, s** - average value of image in square aperture [*x,y*] +-*s* pixels.

**function image\_structure, image, psf** - Returns structure map from *image* and *psf* by formula *image*/(*image* ⊗ *psf*) x *psf*^T, "⊗”=convolution, "^T" = transpose.

**function image\_interp\_flags, image, flags, xi1, xi2, yi1, yi2, dmax** - Interpolate *image* over flagged spaxels, *flags* - 2D data with same x/y axes size as image, with value=1 to be interpolated, value=0 - good pixels, [*xi1:xi2, yi1:yi2*] - co-ordinate range to interpolate over. If not input, then do all spaxels. *dmax* - maximum pixel distance for interpolation (=0 don't test)

**function image\_cut, image, x, y, a** - does **twodcut** at [*x,y*] angle *a* and reset WCS correctly.

**function image\_coord\_map, image** - returns WCS coordinates of corner and central pixels as a 2 x 5 array, with correct projection, using "worldpos" function.

**function image\_trim, image, x, y, w** - returns square image [*x-w:x+w,y-w:y+w*]. Useful for cutouts from a large image.

**function image\_combine, inbuff, mode, omit** - Average/median of multiple images (up to 10). *inbuff* is a multi-extension fits file, created by DPUser _list_(buffer1, buffer2,...). *mode* = 0, average combine, *mode* = 1, median combine. *omit* is value to ignore in DPUser _average_ or _median_ functions. 

**function image\_make\_wcs, wcs_x, wcs_y** - Make a 2D fits file from x and y *wcs* data. wcs data in form [min, max, n]

**function image\_make\_grid, xmin, xmax, xsteps, ymin, ymax, ysteps**  - Make a 2D fits file from *min*, *max* and *steps* in x and y directions. Used to create template for KDE.

**function image\_kde\_gauss, data, grid, kernel** - Create gaussian KDE of x/y data onto 2D map. *data* is 2xn [x,y]. *grid* is 2D image with wcs values set - usually created by **image\_make\_grid** or **image\_make\_wcs**. *kernel* is kernel size as a multiplier of the standard deviations of the data in. If not given, or = 0, then set to Scott's rule, i.e. n^(1/6). If = -1, use Silverman's rule, else use kernel size

**function image_curve_of_growth, inbuff, x, y, r, nrm** - Create curve of growth from *image* for image photometry, centred on pixel *x,y* up to radius *r*. If *nrm* > 0 then normalize result.

### lib\_spectrum
**<u>*Spectrum functions*</u>**

**function spectrum\_make\_disp, val, delt, pix, n** - make 1D vector over range defined by WCS *val*, *delt*, *pix*, *n*.

**function spectrum\_make\_disp\_data, data, axis** - make 1D vector over range defined by WCS values from *data* axis (1,2 or 3)

**function spectrum\_make\_disp\_n, val1, val2, n** - make 1D vector over range [*val1:val2*], number of points *n*

**function spectrum\_mask, spectrum, mask, value, wlflag** - *spectrum* set to *value* between pixel pairs in *mask*. Works for 1D or 3D, assuming last axis is spectrum. *wlflag* =0, mask is in pixels, =1, mask is in wavelength. *value* is usually 0/0 to blank masked sections for polynomial fitting.

**function spectrum\_cont\_slope, spectrum, mask, wlflag** - continuum slope of *spectrum*, masked by wavelength *mask*/*wlflag* (as for **spectrum\_mask**)

**function spectrum\_deslope, spectrum, mask, wlflag** - deslope spectrum, using **spectrum\_cont\_slope** and *mask/wlflag* parameters

**function spectrum\_polyfit, spectrum, order, mask, wlflag** - fit polynomial of *order* to  masked *spectrum* with *mask*, *wlflag*. Returns n x 5 array, 
1 - original data with mask applied
2 - polynomial fit
3 - residual
4 - spectrum-polynomial (continuum subtracted)
5 - spectrum/polynomial (continuum normalised)

**function spectrum\_symm\_flip, spectrum, lambda, part** - split *spectrum* at wavelength *lambda*, flip and add, taking left (*part*=0) or right (*part*=1) sections

**function spectrum\_wave\_to\_lambda, spectrum, l1, l2, nl** - converts a wavenumber *spectrum* to a wavelength spectrum. *l1*..*l2* are a wavelength range to interpolate over with *nl* points. By default, the wavelength range and number of points of the original spectrum are used. WCS values are set.

**function spectrum\_wave\_to\_lambda, wndata** - convert wavenumber spectrum *wndata* to wavelength (nm) with same axis length

**function spectrum\_make\_gauss, spectrum, bi, bs, h, l, w** - make spectrum with gaussian from *spectrum* WCS. *bi, bs* - base intercept and slope, *h* -  height, *lc* - center wavelength, *w* - FWHM (creates artificial gaussian emission line).

**function spectrum\_make\_lorentz, spectrum, bi, bs, h, lc, w** - as for **spectrum\_make\_gauss** but makes a Lorentzian emission line.

**function spectrum\_redisp\_lin, spectrum, data, daxis, xmin,delt, npix, zero, norms,  prnt, fluxcons** - re-disperse a *spectrum*. Parameters are:

- *data* - data with dispersal solution (if =0 then use parameters for dispersion)
- *daxis* - spectral axis of data (default is last axis of *data*)
- *xmin, delt, npix* - dispersion solution if data=0
- *zero* - if =1, then set redispersed spectra to zero where out of original range, rather than NaN (=0)
- *norms* - if =1, normalize dispersed spectra \[0,1\] (default 0)
- *prnt*  - if =1, print spectral range information (default 0)
- *fluxcons*  - if = 1, spectrum is flux, rather than flux density, so conserve total (default 0)

**function spectrum\_from\_xy, spectrum** - re-disperse *spectrum* from 2D x and y bintable to wavelength range and same number of points.

**function spectrum\_from\_tablexy, data, l1, l2, npix, xscl, yscl** - re-disperse *spectrum* from 2D x and y bintable to wavelength range *l1..l2* and number of points *npix*. x and y values are scaled by *xscl* and *yscl* respectively.

**function spectrum\_redisp\_xy, spectrum, l1, l2, npix, xscl, yscl** - as for **spectrum\_from\_tablexy** but from standard *spectrum*.

**function spectrum\_redisp\_data, data, l1, l2, npix, xscl, yscl** - as for **spectrum\_from\_tablexy**, but from image *data*, wavelength in col 1, flux in col2

**function spectrum\_from\_dataxy, xdata, ydata, l1, l2, npix, xscl, yscl** - re-disperse spectrum from 2D x and y data to wavelength range [*w1, w2*] with step *delt*.

**function spectrum\_interp, spectrum, x1, x2** - Smooth over bad pixels [*x1:x2*].

**function spectrum\_sn, spectrum, window** - Estimate spectrum S/N from itself - not 100% accurate but good for comparisons, *window* is smoothing and noise estimation window. Returns vector of same length as *spectrum* with S/N estimate, blank where *spectrum* is 0.

**function spectrum\_sn\_wl, spectrum, spectrum, lambda, width** - estimate spectrum S/N from itself at a defined wavelength (*lambda*) /window width (*width*), i.e. *lambda* ± *width*. Returns [*signal, noise, SN*] (average and standard deviation over the window)

**function spectrum\_clean, inbuff, thresh** - clean *spectrum* using **dpixcreate**/**dpixapply**. If *thresh* is not sepecifed the threshold for **dpixcreate** is set to median(*spectrum*)/2.

**function spectrum\_apply\_snr, spectrum, snr** - as for **[cube\_apply\_snr](#cube\_apply\_snr)**. *snr* can be a single value or spectrum.

**function spectrum\_wave\_to\_vel, inbuff, clambda** - returns *inbuff* with wavelength axis changed to velocity, with zero value at *clambda*.

**function spectrum\_comb\_sigma, data, sigma, toler, omit** - spectrum combine with sigma clipping algorithm. *data* is a 2d array with each row a spectrum. *sigma* is standard deviations allowed (default 3). *toler* is tolerance for convergence on clipping algorithm (default 0.1). *omit* is value to omit on averaging (if parameter not entered, don't omit any value).

**function spectrum\_flux\_ul, inbuff, lambda, width, deslope, ignore, flag** - see **[cube\_flux\_ul](#cube\_flux\_ul)** function description.

**function spectrum_check_monotonic, wbuff, dbuff** - Check spectrum wavelengths in *warray* are monotonic, and remove incorrect ones for *darray*. If *darray* is zero or nor given *warray* has x and y values in cols 1 and 2. Returns is 2 x n array

### lib\_io
**<u>*Input/output to and from text and fits files*</u>**

**function io\_text\_FITS\_1D, bintable2d** - converts string array buffer *bintable2d* with format of "wavelength, data" to spectrum fits data, setting WCS values. Assumes wavelength is evenly spaced. Note **import** function of QFitsView does very similar (with more parameters).

**function io\_text\_FITS\_3D, bintable2d, nx, ny, nz, blank**  - converts string array buffer *bintable2d*, with format of “i,j,v1,v2..." to fits data cube size [*nx,ny,nz*]. Default value for resulting cube is *blank* (e.g. 0 or 0/0) - can have missing [*i,j*].

**function io\_text\_FITS\_interp, fname, xstart, xdelta, xnum, xscale, yscale, ignore** - converts text from file *fname*, with format of "wavelength, data" to spectrum fits data, setting WCS values. The values are interpolated to the range defined by *xstart*, *xdelta* and *xnum*. Wavelength and data value are scaled by *xscale*, *yscale* (default 1). *Ignore* lines at the start are skipped (e.g. column headers).

<a name="io\_FITS\_text\_1D"></a>**function io\_FITS\_text\_1D, spectrum, prefix, cutoff** - converts *spectrum* to text, CSV format, line 1 = “*prefix*\_\_Wavelength, *prefix*\_Counts”. Values below *cutoff* (non-zero) are set to NaN (Be aware of QFitsView Edit > Copy functionality)

**function io\_FITS\_text\_2D, image, prefix** - converts *image* to text in CSV format, line 1 = “*prefix*\_Wavelength, *prefix*\_Flux\_1, *prefix*\_Flux\_2 …. "

**procedure io\_FITS2TXT\_1D, fname, cutoff** - converts 1D FITS to text file *fname* assuming file is in  working directory - output is same as input file with “.txt” type. *Cutoff* as for **[io\_FITS\_text\_1D](#io\_FITS\_text\_1D)**.

**procedure io\_FITS2TXT\_2D, fname** - as above but for image (2D) file

**function io\_cube\_from\_xyz, cube,bintable2d, n** - make a cube from *bintable2d*, *cube* is template, resized to *n* on axis 3,  first 2 values in data are x,y co-ords, rest are values along z axis

**function io\_import\_TXT\_1D, fname** - import data from file *fname* in text format

**function io_col_to_csv, inbuff, sep** - Create CSV line from 1D fits array *inbuff*. *sep* is the element delimiter, default ",". If = "T", use \<TAB\> character. If = "C" or "", use "," character.

**function io\_bintable\_text, fname, extn, sep** - Export binary table to CSV. *fname* - file name of fits bintable. *extn* - extension name or number (default 1). *sep* - see **io\_col\_to\_csv** function. Returns a string array of all lines 

**function io\_writefile\_ctr, inbuff, prefix, suffix** - Write file with *prefix*, *suffix* and counter (limit 9999), without overwriting existing files. For example, 

`io_writefits_ctr data, "TEST_", ".fits" `

will output TEST\_0001.fits, TEST\_0002.fits ... etc. to working directory.


### lib\_masking
**<u>*Masking functions for images and cubes*</u>**

**function mask\_from\_image, image, level, low** - create a mask from data *image*, setting to 1 if > *level*, to *low* (usually 0) if \<*level*

**function mask\_from\_image\_nan, image, zero** - create a mask from data *image*, setting to 1 if  data value<>NaN. If *zero* = 0 or 1, set mask to NaN or 0 at NaN values.

**function mask\_data, image, level, low** - masks data *image*, setting to *low* if < *level*

**function mask\_data\_median, image, level, low** - as above, but sets data image > *level* to median of *image*

**function mask\_circle, data, x, y, r, v, rev** - masks *data* (image or cube) with circle center [*x,y*] radius *r*, set masked-out value to *v* (default 0). If *rev*<>0, reverse mask.

**function mask\_set\_nan\_min, data, minvalue** - set *data* values to *minvalue* if value = NaN. If minvalue is zero, use the current minimum value. Equivalent to **cblank** function is *minvalue*=0

**function mask\_cone, data, xc1, yc1, xc2, yc2, pa, beta, maskflag** - mask cone area over *data* (either image or cube), with equator [*xc1, yc1*], [*xc2, yc2*] (can be same coordinates for a point apex), centerline angle *pa* (from positive x-axis), internal full-angle *beta*. If *maskflag*=0, return the mask, if *maskflag*=1, return the masked input data.

**function mask\_line, data, x1, y1, x2, y2, side** - creates a mask of *data* dimensions on one side of a line [*x1, y1*], [*x2, y2*]. *side* =0 for left, =1 for right side of line

### lib\_velmap
**<u>*Velocity map (velmap) extension functions*</u>**

<a name="velmap\_std\_to\_ext"></a>**function velmap\_std\_to\_ext, velmapstd, r, cmin, vmethod, vzero, vx, vy** - convert standard QFitsView velmap *velmapstd* to extended form, *r*=instrumental resolution, *cmin*= minimum continuum value. Output is in extended velmap format - see below. Velocity zero is set by *vmethod* =

- 0 - median   
- 1 - average
- 2 - flux-weighted average
- 3 - manual (*vzero* value)
- 4 - pixel ([*vx,vy*] is set to zero)

**function velmap\_vel\_center, velmapstd, vmethod, vcenter, vx, vy** - returns the wavelength value from the standard *velmapstd* cube, using the methods as above

**function velmap\_vel, velmapstd, vmethod, vcenter, vx, vy** - returns the velocity map from the standard *velmapstd* cube, using the methods as above

**function velmap\_vel\_set, velmapext, vmethod, vcenter, vx, vy** - Fix extended *velmap* velocity as per **velmap\_std\_to\_ext** (re-do extended velmap cube)

**function velmap\_rescale, velmapext, scale** - rescales extended *velmap* flux data (e.g. flux calib change)

<a name="velmap\_fix"></a>**function velmap\_fix, velmap, contlo, conthi, flo, fhi, vlo, vhi, wlo, whi, snmin, setvalue** - clean up *velmap* (either standard or extended form), setting values out of range to *setvalue*. Value ranges 

- *contlo, conthi* - continuum
- *flo, fhi* - flux
- *vlo, vhi* - wavelength
- *wlo, whi* - fwhm
- *snmin* - signal/noise minumum, detemined by velmap height and height error at each pixel.
- *setvalue* - value to set where spaxel is out of range (default 0/0)

**function velmap\_clean\_pixels, velmap, n** - reverse of velmap\_fix\_interp, removes lone pixels that have NaN as neighbours. *n* is the maximum number of non-NaN neighbours, by default = 0 (i.e. pixel must be surrounded bu NaNs)

**function velmap\_extcorr, velmap, av, lambda** - extinction correct velocity map *velmap* at wavelength *lambda* (in nm), *av*=extinction A\_V

**function velmap\_extcorr\_map, velmap, extmap, lambda** - as above, but *extmap* is a map of extinction values

<a name="velmap\_fix\_interp"></a>**function velmap\_fix\_interp, velmap, npix** - interpolate velmap *velmap* missing values, indicated by NaN in continuum layer (i.e. layer 1) (usually after **[velmap\_fix](#velmap\_fix)**). *npix* is interpolation width maximum

**function velmap\_clean\_map\_wvt, velmap, map, nregion** - Clean up velmap *velmap* based on WVT *map* region number, setting region *nregion* pixels to NaN

**function velmap\_mask, velmap** - set *velmap* to NaN where continuum=0

**procedure velmap\_comps, velmapext, prefix, hmax** - Output *velmapext* components from *velmap*, to the current working directory. *prefix* (string) sets file names, terminated with 

- \_Flux - flux (layer 12)

- \_Flux\_Norm - normalized flux (range [0..1])
- \_Vel - velocity (layer 10)
- \_EW - equivalent width (layer 13)
- \_Sig - dispersion (layer 11)
- \_VelHist - velocity histogram. If *hmax* > 0, a 50 bin histogram of the velocity with range *hmax* -> *hmax*
- \_SigHist - - dispersion histogram. If *hmax* > 0, a 50 bin histogram of the velocity with range 0->*hmax*

**function velmap\_from\_profit, profit\_data** - convert *profit\_data* PROFIT cube format (see Riffel, R. A. 2010, Astrophys Space Sci, 327, 239, http://arxiv.org/abs/1002.1585) for standard velmap format.

**function velmap\_derotate, cube, velmodel, lambdac** - Subtract velocity model *velmodel* from a data *cube*, where the velocity is determined from central wavelength *lambdac*. This shifts each spaxel spectrum by a wavelength amount calculated by the central wavelength and velocity model. This might be used if you have, say, a stellar rotation model and you want to apply it to gas emission lines.

**function velmap\_flux\_fix, velmap, scale** - rescale velmap fluxes (both standard and extended types) by scale. This corrects for e.g. incorrect flux calibration.

### lib\_chmap 
**<u>*Channel map functions*</u>**

<a name="chmap\_create"></a>**function chmap\_create, cube, lambda\_cent, lambda\_width, cutoff, width\_factor, smooth** - make a channel map from the cube

- *lambda\_cent* - estimate of central wavelength
- *lambda\_width* - estimate of FWHM
- *threshold* - % of maximum for cutoff - default 0
- *width\_factor* - wavelength widow (multiple of *lambda\_width*) - default 2.5
- *smooth* - integer=boxcar, non-integer=gauss, 0=no smoothing - default 0

Returns a cube of channel maps, with axis 3 in velocity difference (km/s) from median. Spaxel values are FLUX (not flux density) in that channel.

<a name="chmap\_rebin"></a>**function chmap\_rebin, cube, lnew, velwidth, sm, minval**- rebin channel maps in *cube* into *lnew* bins between velocities *v1* and *v2* (usually symmetric about 0, but not necessarily), with *sm* smoothing value, integer=boxcar, non-integer=gauss, 0=no smoothing, set output to NaN where < *minval*

<a name="chmap\_comps"></a>**procedure chmap\_comps, cube, dirout, fnameout** - splits channel map *cube* into components and writes images to folder *dirout*, named *fnameout* plus velocity (e.g. if *fnameout* ="pa\_beta-450”, then output file name will be e.g.  "pa\_beta-100.fits" etc.

**function chmap\_from\_velmap, velmapstd, cube\_template, width, res** - create a channel map from a standard velmap *velmapstd*. The velmap is evaluated using *cube\_template*, then the channel map is generated using the median velocity from the velmap with width multiplier (how far to extend the channels over the median FWHM) *width* (default 2.5). The FWHM is corrected by spectral resolution *res* (default 0).

### lib\_pv
**<u>*Position Velocity Diagram functions*</u>**

**function pv\_array, cube, ystart, wslit, nslit, lcent, lwidth** - create pv diagram from *cube* parallel to x axis, *ystart* - y pixel to start, *wslit* - slit width, *nslit* - number of slits, extract over range *lcent*-*lwidth* to *lcenter+lwidth*

**function pv\_single, cube,  xc, yc, angle, width, lcent, vwidth, npix, contflag** - extract single PV plot at *xc*/*yc*/*angle*/*width* - centerered on *lcent*. *vwidth* - velocity width around *lcent*, rebinned in velocity to *npix* channels. *contflag* =1 subtract continuum (flux) =2 divide continuum (effectively the same as equivalent width) =0 don’t remove continuum

**function pv\_ratio, cube,  xc, yc, angle, width, lcent1, lcent2, vwidth, npix** - create PV diagram as above for ratio of 2 lines *lcent1*, *lcent2*

**function pv\_meddev, image** - divide *image* by median along x axis (useful for EW for PV diagrams)

### lib\_wvt
**<u>*Weighted Voronoi Tesselation functions*</u>**

**function wvt\_cube, cube, sn\_target** - make WVT *cube* using noise in each spaxel, to *sn\_target*. Bad pixels where S/N is > 10x brightest pixel S/N

<a name="wvt\_cube\_mask"></a>**function wvt\_cube\_mask, cube, l1, l2, mask, cutoff, sn1, sn2** -  make WVT cube using 2 S/N ratios, inside and outside *mask*. Returns WVT applied to *cube*.

- *l1, l2* - wavelength range to use for signal and noise determination (“quiet” part of spectrum with no emission lines)
- *mask* - if 2D mask, use this. If *mask*=0, use *cutoff* to determine mask
- *cutoff* - percentage of peak maximum for mask level
- *sn1*, *sn2* - S/N ratios for inside/outside mask. If *sn2*=0, just use *sn1* over whole cube

**function wvt\_sn\_mask, cube, l1, l2, mask, cutoff, sn1, sn2** - as for **[wvt\_cube\_mask](#wvt\_cube\_mask)**, except returns WVT image data - layers:

1. Signal
2. Noise
3.  S/N
4.  Mask
5.  Signal binned
6.  Signal bin map
7.  Bin density (1=maximum - smallest bins , 0=minimum - biggest bins)

**function wvt\_build\_from\_map\_cube, cube, wvtmap, prntflag** - make WVT cube from *cube* and *wvtmap*. If *prntflag* =1 print diagnostic every 100 regions

**function wvt\_build\_from\_map\_image, image, wvtmap** - make WVT image from *image* and *wvtmap*.

**function wvt\_velmap, velmap, layer, sn** - make WVT velmap from standard or extended velmap *velmap*, *layer* is either 0=continuum, 1=flux, *sn*=S/N target

**function wvt\_density, wvtmap** - make map of region density, i.e. 1/# of pixels in region. *wvtmap* is WVT with /map flag.

**function wvt\_cube\_to\_specarray, cube, wvtmap, normflag, prntflag** - convert *cube* inbuff to spectrum array, using *wvtmap* regions. If *normflag* = 1, divide each spectrum in array by the first one. If *prntflag* = 1, print running diagnostics

**function wvt\_specarray\_to\_cube, image, wvtmap** - reverse of **wvt\_cube\_to\_specarray**

### lib\_general

***<u>Miscellaneous functions</u>***

**function indexreform, index, xsize, ysize, zsize** - returns 3D co-ords from 1D *index*, given dimensions *xsize, ysize, zsize*. Values returned as array.

**function lognan, data** - set log of *data*, setting zero and NaN values to NaN

**function clipnan, data, low, high** - set values outside range [*low..high*] to NaN

**function axiscentroids, image, axis** - returns centroids of each image row/column, row-*axis*=1, column-*axis*=2; used e.g. for finding centroids of pv diagram

**function histogram\_bin, data, low, high, bin, normflag** - create a histogram from inbuff data (any dimensions), from \_low\_ to \_high\_ values in *bin* bins. Histogram is normalised if *normflag* =1. Output x-axis values are set to range.

**function profile\_export, data, scale1, scale2, scale3, offset** - Export 1D profiles from *data* with up to 3 separate scales, e.g. arcsec, pc, Re plus the pixel scale, \_offset\_=1 offsets by 1/2 a pixel (e.g. for log scale plot) (default 0). Scales default to 1 if not given.

**function butterworth\_filter, order, cutoff, size** - Create a Butterworth filter for order \_order\_ for a square of sides \_size\_, with \_cutoff\_ Nyquist frequency.

**function interp, data, x, x1, x2**  - linearly interpret over *data* at position *x* over *x1*-*x2*.
Used to image\_interp\_flags and cube\_interp\_flags

**function fits\_round, inbuff** - returns *inbuff* with all values rounded.The ESO pipeline for the SOFI instrument required this for raw images - if these had been manually flat-fielded, the pipeline failed.

**function inf\_clean, inbuff** - returns *inbuff* with all values of "1/0" (i.e. Inf) set to 0/0 (i.e. NaN).

### lib\_astronomy

**<u>*General astronomy functions, implemented from IDL.*</u>**

**function G** - gravitational constant (MKS)

**function Msun** - Mass of the sun in kg

**function Pc** - 1 Parsec in meters

**function airtovac, wave** - Convert air wavelengths to vacuum wavelengths, *wave* in Å

**function planck, wave, temp** - Calculates the Planck function in units of ergs/cm2/s/Å. *wave* in Å, *temp* in degrees K.

**function coordstring, ra, dec, rad, format** - create a nice string of celestial coordinates. *ra* and *dec* should be given in radians; if *rad*>0, convert ra and dec to radians (default ra and dec input in degrees). *format* sets the output format, 0 = long string, 1 = short string array [ra, dec], 2 = number array [rah, ram, ras, ded, dem, des].

 Example:

`>>>print coordstring(9.1234,5.678)`=>`RA = 0h 36m 29.62s, DEC = + 5d 40' 40.8"`
`>>>print coordstring(9.1234,5.678,0,1)` =>
`00:36:29.616`
`+05:40:40.80`
`>>>print coordstring(9.1234,5.678,0,2),/values`=>
`0`
`36`
`29.616`
`5`
`40`
`40.8`

**function sextodec, sexstr** - converts sexadecimal string to decimal number, e.g.
`print sextodec("01:02:03")`=>` 1.03417`

**function dectosex, decnum** - reverse of **sextodec**, e.g. 
`print dectosex(1.2345)`=>`1:14:04.2`

### lib\_astro\_general

**<u>*General astrophysics functions.*</u>** 

All the "lib\_astro\_*.DPUser" functions are executed from the "lib\_astro.DPUser" script.

**function redshift\_data, data, z, rdflag, nanflag, smth, fcons** - redshift *data* by *z*, assuming last axis is wavelength (this works on spectra or cubes). WCS values are set. If *rdflag* = 1, redisperse shifted data to same wavelength range as input. If *nanflag*=1, set NaNs in output at same pixels as in input. If *smth*>0, gaussian kernel smooth by no. of pixels. If *fcons*=1, conserve total values.

**function bb, temp, wl** - black-body value for temperature *temp* and wavelength *wl* (in m)

**function bb\_make,temp, wl1, wl2, npix, wlflag** - make black-body function at temparture *temp*, wavelength range *l1* to *l2*, number of pixels *npix*, *wlflag* sets the wavelength scale of *wl1* and *wl2*.
0 = Å
1 = nm (default)
2 = μm

**function bb\_make\_log, temp, wl1, wl2, npix, wlflag** - make bb at temp *temp* over log wavelength [*wl1,wl2*], creating spectrum length *npix*, wavelength scale flag *wlflag*

**function bb\_div, spectrum, temp, wlflag** - divide *spectrum* by black-body at temperature *temp* and wavelength scale flag *wlflag*.

**function diag\_ne\_SII, lSII1, lSII2, Te** - Electron density diagnosis from [SII] line ratio. *lSII1/lSII2* are flux/luminosity of [S II] 6717/6730Å emission lines, *Te* is electron temp. defaults to 10000K. Formulae from _Proxauf+ 2014_.

**function diag_Te_OIII, lOIII1, lOIII2, lOIII3** - electron temperature from *lOIII1,2,3* are flux/luminosity of [O III] λ 5007/ 4958/4363Å. Formulae from _Proxauf+ 2014_.

### lib\_astro\_mapping

**<u>*Astronomy functions (mapping and excitation diagrams)*</u>**

**function map\_compare\_diagram, image1, image2, min1, max1, min2, max2, nbin, lgaxesflag** - Map diagram density plot. *image1*, *image2* - value maps, x and y axes. *min1/2*, *max1/2* - min and maximum values for axes 1/2. *nbin* - no of bins on each axis. *lgaxesflag* - 1=plot in log space (min,max must be in log values). Generates e.g. BPT diagrams

**function map\_compare\_pos, image1, image2, image3, image4, x, y, boxsize** - get 2 sets of map ratios (*image1*/*image2*, *image3*/*image4*)at position [*x,y*], averaged over *boxsize* x *boxsize* pixels (e.g. excitation ratios at feature position)

**function map\_basis\_distance, basex0, basey0, basex100, basey100, x1, x2, y1, y2, size** - creates an image (dimensions *size*, limits (*x1, y1*), (*x2, y2*)) of distance from basis points 0 to 100% [*basex0, basey0*] to [*basex100, basey100*] - for use in AGN mixing ratios for contour values.

**function map\_compare\_basis, image1, image2, basex0, basey0, basex100, basey100, lgaxesflag** - plots basis distance (AGN mixing ratio) from basis points [*basex0, basey0*] to [*basex100, basey100*] . *lgaxesflag* - 1=take log of *image1*, *image2* before calculation

**function map\_regime\_ir, image1, image2, a1, a2, a3, b1, b2** - create position excitation map. *a1, a2, a3* are the SF, AGN, LINER, b1, b2If a1=0, use the standard infrared _Riffel 2013_ excitation regimes. *image1* is H$_2$/Brγ, *image2* is [Fe II]/Paβ. Both images are in log values. Output values at each spaxel are SF=1, AGN=2, LINER=3, TO1=4, TO2=5

**function map\_regime\_optical, image1, image2, typeflag**- create position excitation map for optical line ratios (*image1* and *image2*) from _Kewley et al. 2006_ regimes. For *image1*, *typeflag* = 1 ([N II]/Hα diagram), =2 ([S II]/Hα diagram), =3 ([O I]/Hα diagram). *image2* always has [OIII]/Hβ. Returns at each pixel 1=SF, 2=Seyfert, 3=LINER, 4=Composite.

### lib\_astro\_spectrum
**<u>*Astronomy functions (spectrum)*</u>**

**function spec\_fluxdens, spectrum, l1, l2, prflag** - flux density (counts/nm) for *spectrum* between *l1* and *l2* <u>wavelength</u>; returns a single value. If *prflag*<>0, print results as well.

**function spec\_sn, spectrum, l1, l2** - Compute S/N for *spectrum* over wavelength region *l1*-*l2*, using median rather than average, as more robust.

**function spec\_wave\_to\_vel, spectrum, lambda -** Convert *spectrum* wavelength axis to velocity, with zero velocity wavelength *lambda* .

### lib\_astro\_image
**<u>*Astronomy functions (image)*</u>** 

**function img\_aphot\_annular, image, xcen, ycen, r, ib, ob** - aperture photometry on *image*,centered on [*xcen,ycen*]; aperture *r*, background annulus from *ib* to *ob* (inner to outer boundary). If *ib* and *ob* are zero, set to *r* and 2*r* by default.

**function img\_apphot\_simple, image, xcen, ycen, r, pixsize, scale** - simple aperture photometry on image,centered on  [*xcen,ycen*]; aperture *r*.

**function img\_flux\_to\_mag, image, zpm, ssize, zpflag** - convert flux *image*  to mag, *zpm* is either zero-point magnitude (*zpflag*=0) or zero-magnitude flux (*zpflag*=1 - default). *ssize* = pixel size in arcsec to convert to mag/arcsec^2 (default 1). This can also work for a single number, in which case *ssize* is set to 1.

**function img\_convert\_filters, image1, image2, coeffs12, coeffs21, tol** - Convert images in 2 filters to another filter set using the methodology of ﻿_Holtzman, J. A., Burrows, C. J., Casertano, S., et al. 1995, PASP, 107, 1065_. The coefficients (*coeff12* and *coeffs21*) for WFPC2/WFC3 are from Holtzmann, for ACS from﻿ _Sirianni, et al. 2005, PASP, 117, 1049_. The tolerance for the result is where the maximum magnitude change on each iteration is less than *tol*. 

Example - for conversion from ACS (WFC) F606W and F814W images to *V−I* colour image, the coefficient sets (Table 22 of _Sirianni et al._) are:
*coeffs12* - [26.325, 0.236, 0.000]
*coeffs21* - [25.485, -0.002, 0.000]

The function returns a cube with 3 layers, (1) filter result 1 e.g. *V* (2) filter result 2 e.g *I* (3) difference e.g. *V-I*, as well as printing iteration and result diagnostics.

**function img\_SFD\_dust\_pos, dustimgn, dustimgs, lat, long** - returns galactic dust extinction value, as per _Schlegel et al. 1998, ApJ, 500, 525_. *lat* and *long* are galactic co-ordinates. *dustimgn* and *dustimgs* are the dust mappings, north and south. The standard for *E(B-V)* is the "SFD\_dust\_4096" maps, but others can be used - these are all downloadable from https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/EWCNL5.

**function img\_SF\_dust\_pos, dustimgn, dustimgs, lat, long** - as above but corrects the *E(B-V)* as per _Schlafly, E. F., & Finkbeiner, D. P. 2011, Astrophys J, 737, 103_.

**function img_mag_2MASS, inbuff** - convert 2MASS images to mag/arcsec^2 using img_flux_to_mag and FITS header ZPT. This assumes 1" pixels.


### lib\_astro\_cube
**<u>*Astronomy functions (cube)*</u>** 

**function cube\_apspec, cube, ox, oy, or, bx, by, br, br2, mask** - Get star spectrum from *cube* withusing circular aperture and background, plus a mask circle. 

- *ox, oy, or* - center and radius of aperture
- *bx, by, br* - center and radius of mask, if not required then bx = 0
- *br2* - inner radius of mask annulus (default 0). This annulus allows the background to be centered the same as the aperture (in which case *br2* > *or*)
- *mask* - any other mask required (2D fits) (default none)

Output is spectrum of aperture less average of background (with mask) - values < 0 are set to NaN.

**function cube\_fluxdens, cube, l1, l2, prflag/function spec\_fluxdens, spectrum, l1, l2, prflag** - flux density (counts/nm) between *l1* and *l2* wavelength; returns image (cube\_) or single (spec\_) value. If *prflag*<>0, print results as well

**function cube\_sky\_rem, cube, bckgnd\_lvl** - removes skylines from *cube*. Takes background pixels as those with median value below *bkgnd\_lvl*

**function cube\_sl\_clean, cube, skyline\_list, width** - removes skylines from *cube* using skyline\_linelist (array of wavelengths) - interpolated over wavelength ± *width*

<a name="cube\_clean\_bp\_fix"></a>**function cube\_clean\_bp\_fix, cube, bp\_cube** - cleans *cube* based on bad pixel cube *bp\_cube* using **dpixapply** over x image slices

**function cube\_clean\_bp, cube, threshold** - create bad pixel cube using *threshold* scanning over wavelength slices.

**function cube\_clean\_bp\_limits, cube, ll, ul** - create bad pixel cube from *cube* for input to **[clean\_cube\_bp\_fix](#clean\_cube\_bp\_fix)**, flagging pixels below *ll* and above *ul* values

### lib\_astro\_kinematics

<u>***Astronomy functions - kinematic models***</u>

**function velmodel\_plummer, x0, y0, M0, Re, psi0, inc, xsize, ysize, pixscale, angscale** - create a Plummer kinematic rotation model velocity field in km/s. This can be applied to a cube by e.g. **[cube\_velocity\_correct](#cube\_velocity\_correct)**.
*x0, y0* - centre pixel position
*M0* - enclosed mass (units of Msun)
*Re* - length scale (pc)
*psi0* - line of nodes (default 0 - counterclockwise from +ve X axis)
*i* - inclination of disk (0-90, 0 = face-on - default)
*xsize, ysize* - image size (by default, twice x0, y0)
*pixscale* - pixel scale in arcsec (default 1)
*angscale* - angular scale in pc/arcsec (default 1)

**function velmodel\_hernquist, x0, y0, M0, Re, psi0, inc, xsize, ysize, pixscale, angscale** - as for **velmodel\_plummer** for Hernquist kinematic rotation model.

**function velmodel\_geom, x0, y0, psi0, inc, xsize, ysize** - computes the geometric components of rotational models, with parameters as above. Returns data [*xsize, ysize*, 2] - first layer is rotational component, second layer is radial component. These are independant of the particular model form.

### lib\_astro\_extinction

**<u>*Astronomy functions - extinction*</u>**

**function astro\_extn\_calc, f1, f2, l1, l2, rat, galext, s, flmin1, flmin2**- create an extinction map from 2 emission line maps - *f1*, *f2*  are flux maps, *l1*, *l2*=wavelengths, *rat*=expected flux ratio, *galext*=galactic extinction, *smth*= smooth pixels, *flmin1*, *flmin2*=minimum flux value for each map - calculates the extinction constant (CCM laws for IR and optical)

**function astro\_extn\_correct, data, ebv, law, lscale** - correct *data* for extinction (*ebv* = single extinction value). *data* can be a spectrum or cube, with wavelength from either axis 1 or 3. *law* = 0 -> CCM, = 1 ->Calzetti+00. *lscale* (default 1) - multiplier to convert spectrum wavelength to nm e.g. A->nm = 0.1

**function astro\_extn\_correct\_cube, data, lambda, ebv, law, lscale** - correct value/image *data* for extinction *ebv* at wavelength *lambda*: can be used on single value or an image.

**function astro\_extn\_correct\_spectrum, spectrum, ebv, law, lscale** - Correct *spectrum* for extinction using single  E(B-V)  (*ebv*) value.  

**function astro\_extn\_correct\_cube, cube, ebv, law, lscale** - Correct *cube* for extinction using single  E(B-V)  (*ebv*) value.  

**function astro\_extn\_correct\_lambda, data, lambda, ebv, law** - Correct *data* for extinction at a single wavelength and E~B-V~ (*ebv*) value.

**function astro\_extn\_al\_ccm, lambda, rv** - Calculate _Cardelli, Clayton, and Mathis (1989 ApJ. 345, 245)_ "CCM" extinction curve (A<sub>λ</sub>) at wavelength *lambda* (nm) with R<sub>V</sub> value of *rv* (default 3.1). This includes the update for the near-UV given by _O'Donnell (1994, ApJ, 422, 158)_ 

**function astro\_extn\_al\_cal, lambda, rv** - As for **astro\_extn\_al\_ccm**, but for _Calzetti etal. (2000, ApJ, 533, 682)_ "CAL" extinction curve; *rv* defaults to 4.05.

**function astro\_extn\_ccm, lambda, ebv, rv** - Calculate actual extinction factor for CCM at wavelength *lambda* (nm), E(B-V) *ebv* (default 1) and R<sub>V</sub> *rv* (default 3.1).

**function astro\_extn\_cal, lambda, ebv, rv** - As for **astro\_extn\_ccm**, for "CAL". *rv* defaults to 4.05.

**function function astro\_extn\_const\_ccm, lambda1, lambda2** - Calculate constant in extinction law between two wavelengths (*lambda1, lambda2*) (nm) for  E(B-V) "CCM" calculation.

**function astro\_extn\_const\_cal, lambda1, lambda2** - as for function **astro\_extn\_const\_ccm** for "CAL".

### lib_templates

**<u>*File processing templates*</u>**

**procedure procfiles** - scan and process all files in a directory. 

**procedure procparams** - process files based on a parameter file (CSV).

# Standard Procedures

## Velocity maps
* Create QFitsView **velmap** with wavelength, FWHM estimates
* Examine velmap for continuum, height, wavelength and fwhm “sensible” ranges
* Use **[velmap\_fix](#velmap\_fix)** to clean up velmap, entering "good" ranges for each fit component. This sets out-of-range spaxels to NaN.
* Use **[velmap\_fix\_interp](#velmap\_fix\_interp)** to interpolate over NaN values (if required)
* Use **[velmap\_std\_to\_ext](#velmap\_std\_to\_ext)** to create extended velmap format

**<u>*Standard VELMAP Format*</u>**

1. Continuum
2. Peak height above continuum
3. Wavelength
4. FWHM
5. e\_Continuum
6. e\_Peak
7. e\_Wavelength
8. e\_FWHM
9. Chi-squared

**<u>*Extended VELMAP Format*</u>**

1. Continuum
2. Peak height above continuum
3. Wavelength
4. FWHM
5. e\_Continuum
6. e\_Peak
7. e\_Wavelength
8. e\_FWHM
9. Chi-squared
10. Velocity (zero-point calculated/set by *method* parameter in the **[velmap\_std\_to\_ext](#velmap\_std\_to\_ext)** function)
11. Dispersion (sigma) velocity, corrected for spectral resolution
12. Flux (Peak\*FWHM*1.0699)
13. Equivalent width (flux/continuum)
14. Total support (order + turbulence) (√)
15. Order vs turbulence (|V/σ|)

## Channel maps

- Create basic channel map using **[chmap\_create](#chmap\_create)** (usually do not smooth)
- Rebin to required # of channels (e.g. 9 or 16) using **[chmap\_rebin](#chmap\_rebin)** (smoothing if required)
- Output individual channel maps using **[chmap\_comps](#chmap\_comps)**
<<<<<<< HEAD
# 
=======
>>>>>>> parent of 9e7707d (Minor)
