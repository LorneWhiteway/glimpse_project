# glimpse_project

## Introduction

This repository contains code to allow [glimpse](https://github.com/CosmoStat/Glimpse) (a mass-mapping algorithm that works on patches of the sky that are assumed to be small and flat) to be run on a larger patch (such as the DES footprint). For more information about glimpse see the [glimpse paper](https://arxiv.org/abs/1603.01599).

The code does this in three steps:
1. *create_cutouts*: an input weak-lensing catalogue is subdivided into many smaller overlapping sub-catalogues (call each one a 'cutout');
2. *run_glimpse*: glimpse is run (possibly in parallel) on each cutout;
3. *merge*: the separate glimpse results are merged together to create an output convergence map (in healpix format).

The code is structured so that other mass-mapping algorithms besides glimpse can be handled in the future.

## Interface

The primary interface is the Python script `glimpse_on_curved_sky.py` in the `./scripts` directory. Run this with command line arguments as follows:

| Option | Meaning |
| --- | --- |
| -i | Name of configuration file. See below for details of configuration file format. Required. |
| -t | Task. One of "create_cutouts", "run_glimpse" or "merge". Required. |
| -j | Job control. See below for more information. Optional; see below for default behaviour if not provided. |

Example:
```
./scripts/glimpse_on_curved_sky.py -i ./output/my_job.ini -t create_cutouts
```



### Configuration file

The configuration file should be in [standard INI file format](https://en.wikipedia.org/wiki/INI_file). Below we describe the various sections that should be provided.

#### Section [project]

| Section | Key | Value |
| --- | --- | --- |
| project | glimpse_executable | Path to glimpse executable e.g. /home/Glimpse/glimpse. Required. |

#### Section [create_cutouts]

The *input_catalogue* should be a weak-lensing catalogue. It should be in FITS format, and have (at least) data for RA (decimal degrees [0,360]), DEC (decimal degrees [-90,90]), two shear components, and redshift. Note that redshift is required but is not actually used; if *input_catalogue* does not have redshifts then augment this catalogue by appending a redshift field with dummy values.

Cutout catalogues are FITS files generated by extracting data from *input_catalogue*; use *ra_name*, *dec_name*, *shear_names* and *other_field_names* to specify which columns will be copied from the input catalogue to the cutouts. The field names in the cutout catalogues will be the same as the field names in the input catalogue (i.e. we do not standardise field names).

Each cutout catalogue *c* will contain data on galaxies that are in a diamond-shaped region on the celestial sphere that:
1. is centered on a point that is the centre of a healpixel *h*. Use *nside* to specify the NSIDE for the healpixelisation that sets these centres. The cutout *c* will be assigned an id that is the healpixel id (in RING orientation) of *h*.
2. has side length *cutout_side_in_degrees*;
3. is rotated 45 degrees with respect to the RA/DEC axes.

Once these galaxies have been selected, all the data is then:
1. translated to be centred at the standard centre (RA=180, DEC=0). (We use a standard centre so that we do not need separate glimpse ini files for each cutout; this particular standard centre was chosen as it is well away from the RA=0 line which glimpse does not handle well (due to the discontinuity in RA));
2. rotated (TODO: which direction?) by 45 degrees (so that its sides are parallel to the RA/DEC coordinate axes). Shear values are transformed appropriately. (This rotation is done because our use of healpix to generate cutout centres leads to a preference for diamond-shaped regions, while glimpse prefers a square input region).

If *c* contains no galaxies then no corresponding file will be created.

| Section | Key | Value |
| --- | --- | --- |
| create_cutouts | input_catalogue | The full name of the source weak-lensing catalogue (FITS format). Required. |
| create_cutouts | ra_name | Catalogue field name for right ascension. The values in this field will be written to the cutout files. Optional; default is "RA". |
| create_cutouts | dec_name | Catalogue field name for declination. The values in this field will be written to the cutout files. Optional; default is "DEC". |
| create_cutouts | shear_names | Comma-delimited list of catalogue field names for weak-lensing shear values to be written to the cutouts. Provide one or more pairs of field names. Required. Example: "E1,E2"|
| create_cutouts | other_field_names | Comma-delimited list of catalogue field names for other fields (such as redshift) to be written to the cutouts. Required. |
| create_cutouts | nside | The cutouts will be centred on the centres of healpixels with this value as NSIDE. Optional; default is 16. |
| create_cutouts | cutout_side_in_degrees | The side length of each cutout, in degrees. Optional; default is 16. |

#### Section [merge]

Each glimpse output file contains convergence (kappa) values on an lattice of points (an evenly-spaced square lattice of points in the tangent plane (i.e. tangent at the centre of the glimpse input region) which are then projected back down to the celestial sphere using an [orthographic](https://en.wikipedia.org/wiki/Orthographic_projection_in_cartography) projection).

Merging begins by translating each output lattice of points from the standard centre back to the appropriate healpix centre, and undoing the 45 degree rotation. For each point *q* in a given lattice we assign a weight *w(q)* as follows:
1. if the distance from *q* to the edge (in units given by the lattice spacing; by 'edge' we mean the outermost set of lattice points) is less than *outer_border* then *w(q)* is zero;
2. if this distance is greater than or equal to *inner_border* then *w(q)* is one;
3. otherwise *w(q)* varies smoothly between these two extremes.

We then create a fine healpixelisation with NSIDE *intermediate_nside*. For each pixel *p* in this healpixelisation and each output glimpse lattice *a* we find the point *q(a,p)* in *a* that is closest to the centre of *p* (together with its associated weight *w(a,p)*). The pixel *p* is then assigned a kappa value that is the weighted average of the glimpse kappa values at the points *q(a,p)* (as *a* ranges across all glimpse output files), using as weights the *w(a,p)*. Note that since the weights are zero on the edge of *a*, lattices that are distant from *p* will not contribute to the weighted average (as expected). This gives us a healpix map of weighted average kappa values (we also create a healpix map of weights).

Two final post-processing steps are performed:
1. The value and weight maps are downsampled to NSIDE *output_nside*;
2. If desired (via *apply_galaxy_mask?*) the downsampled value map (but not the weight map) is masked by setting to zero the kappa value in any healpixel in which no source galaxies were found.

| Section | Key | Value |
| --- | --- | --- |
| merge | outer_border | Used to specify weights to be assigned to glimpse lattice points; see above for details. Must be positive. Typical value is 90. Required. |
| merge | inner_border | Used to specify weights to be assigned to glimpse lattice points; see above for details. Cannot be less than inner_border. Typical value is 110. Required. |
| merge | intermediate_nside | NSIDE of intermediate (fine) healpixelisation when merging; see above for details. Typical value is 2048. Required. |
| merge | output_nside | NSIDE of healpixelisation used in final output map. Typical value is 1024. (Fun fact: on Earth an NSIDE of 1024 would give pixels about two-thirds the size of Manhattan). Required. |
| merge | apply_galaxy_mask? | If True, then mask the output map by setting kappa values to zero for pixels containing no source galaxies. Required. |

#### Section [survey]

This section is used by Glimpse to describe the format of input data (which in our case must equal the format of the cutout catalogues).

| Section | Key | Value |
| --- | --- | --- |
| survey | center_ra | This must be set to 180.0 (the RA of the standard centre as described above). |
| survey | center_dec | This must be set to 0.0 (the DEC of the standard centre as described above). |
| survey | size | The side length in degrees of the glimpse output lattice of points. It is optimal to set this to the same value as create_cutouts::cutout_side_in_degrees (which is typically 16). |
| survey | units | Set this to 'degrees'; if another unit is chosen then adjust values in this section accordingly. |
| survey | hdu | Set this to 1 (this corresponds to the format of the cutout catalogue files). |
| survey | flip_e2 | Set to true or false depending on the weak lensing shear quote convention used in the input catalogue. Experiment to find the correct value. TODO: verify that we handle properly the 'rotate by 45 degrees' calulation for shear when this setting is 'true' (it works OK when 'false'). TODO: Describe the two quote conventions in more detail. |
| survey | ra | Should equal create_cutouts::ra_name. |
| survey | dec | Should equal create_cutouts::dec_name. |
| survey | e1 | The name of the first shear field (which must have been mentioned in create_cutouts::shear_names). |
| survey | e2 | The name of the second shear field (which must have been mentioned in create_cutouts::shear_names). |
| survey | z | The name of the redshift field (which must have been mentioned in create_cutouts::other_field_names). |

#### Section [cosmology]

This section is used by Glimpse to control the fiducial cosmological model. TODO: Suspect that this is not paid attention to as we are not using source redshifts?

| Section | Key | Value |
| --- | --- | --- |
| cosmology | Omega_m | Matter energy parameter (e.g. 0.25) |
| cosmology | h | Hubble parameter in units of 100 km/s/Mpc (e.g. 0.70). |


#### Section [field]

This section is used by Glimpse in part to set the size and scaling of the glimpse output lattice.

| Section | Key | Value |
| --- | --- | --- |
| field | units | Set this to 'arcmin'; if another unit is chosen then adjust values in this section accordingly. |
| field | pixel_size | The grid spacing for the glimpse output lattice, in arcmin. Typical value is 3.5. |
| field | padding | Number of rows of extra 'padding' lattice elements at the edges of the glimpse output lattice. Typical value is 28. |
| field | include_flexion | Set this to false (so that second order 'flexion' information is not used). |
| field | zlens | Set this to -1 (so that source galaxy redshifts are not used). |

With a pixel size of 3.5 arcmin, an output lattice of 16 degrees on a side, and 28 padding elements, we get a glimpse output lattice of 330 points, squared, covering 256 sq deg (so 425 lattice points per sq deg). With outer_border and inner_border set to 90 and 110 respectively, we essentially keep only the central one-ninth of each lattice i.e. 28 sq deg. With the cutout centres based on a healpix NSIDE of 16, each cutout is centred on healpixel of 13 sq deg, so the overlap ratio is 2.1 (= 28/13). Thus the resulting lattice point density is 900 lattice points per sq deg i.e. 0.25 lattice points per sq arcmin. With an output map NSIDE of 1024 each output pixel is 11.8 sq arcmin, for a density of 2.95 lattice points per output pixel.

#### Section [parameters]

This section controls the algorithm used by Glimpse.

| Section | Key | Value |
| --- | --- | --- |
| parameters | nrandom | Unclear; suggest keeping this at the glimpse example value of 1000. |
| parameters | niter | The number of iterations to use in glimpse's 'primal_dual' algorithm. Suggest values of 500 or 1000. TODO: Investigate consequences of non-convergence of this algorithm. |
| parameters | nreweights | Reweighting is the first of two procedures to correct bias in the amplitude of the output values; see section 3.2 in the [glimpse paper](https://arxiv.org/abs/1603.01599). Niall Jeffrey recommends setting this value to 0 to turn off this procedure. |
| parameters | niter_debias | Reweighting is the second of two procedures to correct bias in the amplitude of the output values; see section 3.2 in the [glimpse paper](https://arxiv.org/abs/1603.01599). TODO: This value needs to be explored; to date it has been set to 0. |
| parameters | nscales | Controls the number of wavelet scales used by glimpse. Section 5.2 of the [glimpse paper](https://arxiv.org/abs/1603.01599) describes an example in which this parameter is not crucial. Suggest keeping this at the glimpse example value of 7. TODO: experiment with alternative values. |
| parameters | lambda | Regularisation parameter, controlling the tradeoff between likelihood (i.e. finding a reconstruction that explains the data well) and prior (i.e. finding a reconstruction with high prior probability i.e. is sparse). Niall Jeffrey suggests using 3.0 for observed data (this was found optimal on DES SV data - see section 4 in [this paper](https://arxiv.org/abs/1801.08945)). The [glimpse paper](https://arxiv.org/abs/1603.01599) (section 3.5) suggests using 0.01 for simulated data in which shape noise is not present. |
| parameters | battle_lemarie_reg | Unclear; suggest keeping this at the glimpse example value of 1. |
| parameters | last_scale_reg | Unclear; suggest keeping this at the glimpse example value of 2. |


### Job Control

Use this command line parameter to subselect only some data for processing.

1. When creating cutouts and when merging, use this to specify that only a subset of all possible cutouts should be handled - this is useful primarily during testing. Use Python slice notation, so that for example `[2000:2100]` would handle only cutouts with 2000 <= id < 2100, while `[::2]` would handle all even-numbered cutouts. Note that ids are zero-based. For cutout creation and merging, job_control is optional; if omitted then all possible cutouts will be handled.
2. When running glimpse, use this to specify the id of which single cutout is to be processed; this will be useful for example in a job script used to run glimpse in parallel. Note that ids are zero-based. In glimpse mode this parameter is required.

## File names

1. At run time the user specifies an ini file. All intermediate and output files are created in the same directory as the ini file (so it makes sense to give this directory a name reflecting its purpose).
2. The ini file can have any name (but glimpse.ini is traditional).
3. Cutout catalogues will have names `{id}.cat.fits` where `{id}` is the id of the cutout (= healpix id of central point).
4. Output files from glimpse will have names `{id}.glimpse.out.fits`.
5. Final output value and weight files will be called `glimpse.merged.values.dat` and `glimpse.merged.weights.dat` respectively.


## Information useful for running glimpse on the splinter cluster at UCL.

### To build Glimpse in the splinter environment

1. git clone https://github.com/LorneWhiteway/glimpse_project.git
2. cd ./glimpse_project
3. source ./set_environment
4. git clone https://github.com/CosmoStat/Glimpse.git
5. Edit ./Glimpse/src/field.cpp to fix a bug: Find the line 
	`lensKernel[ind] = 1. ;`
and in the same `for` loop also add the line
	`lensKernelTrue[ind] = 1. ;`
6. cd ./Glimpse
7. ../cm.sh (This will run cmake with the correct settings)
8. cd ./build
9. make

### Things to know about glimpse:
1. Command line arguments should be ini file, input catalogue, and output file name (in that order), The input item specifiers '--config', '--data' and '--output' should not be provided (despite what the function help says).
2. In the ini file, the 'size' argument is the length of the edge of a square. 
