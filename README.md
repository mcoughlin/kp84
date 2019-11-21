# kp84
Kitt Peak 84 inch scripts</br>
Instrument paper: [Coughlin et al. 2019](https://arxiv.org/abs/1901.04625)

## Usage
The following list steps to reduce photometric data, using 2019-11-17 as an example.

### `kp84_setup_reduction.py`
`python kp84_setup_reduction.py --day 20191117`
1. Basic calibration steps: create master bias, dark, and flat frames.<br>
For calibration files done in mode 0, they need to be downsampled by 2x2.<br>
Although master bias file is not used, as it is essentially the same as the master dark in mode 0.<br>
Typically, flats are taken with filter sloan _gr_ and Johnson _(U)BVRI_.
2. Processing science frames (subtract dark, divide flat).
3. Make register folder; Solve astrometry and save the wcs, using [astrometry.net](http://astrometry.net/).<br>
Call `kp84_get_wcs.py`.<br>
The wcs in the original fits header is often off by 1--2 arcmin.<br>
Shifts between each frames in the multi-extension cubes are calculated in this step and saved to the registration folder.
- Default upload image is the best frame in each cube (the one with most point sources identified). Specify central (ra, dec) with uncertianty = 5 degrees.<br>
- If astrometry failed after trying 3 minutes, then stack all images, using the first extension as referencce.<br>
I took the median of un-shifted region, try 3 minutes, specify central (ra, dec) with uncertianty = 10 arcmin this time.
- If astrometry still fails using the stacked image, then the object's position (x, y) must be given to the following script.

### `kp84_photometric_reduction.py`
`python kp84_photometric_reduction.py --day 20191117 --objName ZTFJ01395245 --doDifferential --doSaveImages`
1. Find the coordinate of object (from the file `input/observed.dat`). So make sure to add this beforehead.<br>
Then for each file (`kped_20191117_hhmmss_ZTFJ01395245*_cl_o`) that belong to the object ZTFJ01395245, do the following steps:
2. Use the wcs, find the (x, y) of object in each frame, save to the processing fits file's headers<br>
Mask frames where the object shifted outside of the field.
3. Photometry
- Copy the pre-processed image into output directory, name it as `science.fits`
- Run [Source Extractor](https://www.astromatic.net/software/sextractor) to identify point sources. <br>
`sex science.fits -c default.sex -PARAMETERS_NAME daofind.param -FILTER_NAME default.conv -CHECKIMAGE_TYPE BACKGROUND -CHECKIMAGE_NAME science.background.fits -CATALOG_NAME science.cat -MAG_ZEROPOINT 0.0`</br>
All default files are in the `/defualt` directory. 
- Run forced photometry using [PythonPhot](https://github.com/djones1040/PythonPhot/blob/master/PythonPhot/aper.py)<br>
The default aperture size is 10 pixels. This can be changed by setting `--aper_size`


Some notes:
- If transients, turn on `--doSubtraction --subtractionSource ps1`, then `SExtractor` will also run on `science.sub.fits`.
- If there are enough objects to solve for zero point, turn on `doZP`<br>
This can be hard sometimes due to the limited field of view (4x4 arcmin)
- If too faint, then turn on `--doStack --nimages 5`
- If do not turn on `--doSaveImages`, then the `science.fits` file will be deleted after photometric reduction.
