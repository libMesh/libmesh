Software for Manipulating or Displaying NetCDF Data {#software}
===================================================

[TOC]

This document provides references to software packages that may be used for manipulating or displaying [netCDF](/software/netcdf/) data. We include information about both freely-available and licensed (commercial) software that can be used with netCDF data. We rely on developers to help keep this list up-to-date. If you know of corrections or additions, please [send them to us (mailto:support@unidata.ucar.edu). Where practical, we would like to include WWW links to information about these packages in the HTML version of this document.

Other useful guides to utilities that can handle netCDF data include ARM's list of [ARM-tested netCDF data tools](http://science.arm.gov/%7ecflynn/ARM_Tested_Tools/), which includes some downloadable binaries and the NOAA Geophysical Fluid Dynamics Laboratory [guide to netCDF utilities](http://nomads.gfdl.noaa.gov/sandbox/products/vis/data/netcdf/GFDL_VG_NetCDF_Utils.html).

------------------------------------------------------------------------


Freely Available Software {#freely}
=========================

ANDX and ANAX {#ANDX}
------------------------------------

The ARM Program has developed [ANDX (ARM NetCDF Data eXtract)](http://engineering.arm.gov/~sbeus/andx-web/html/), a command-line utility designed for routine examination and extraction of data from netcdf files. Data can be displayed graphically (line-plot, scatter-plot, overlay, color-intensity, etc.) or extracted as ASCII data. Whether displayed graphically or extracted as ASCII, results can be saved to disk or viewed on screen.

[ANAX (ARM NetCDF ASCII eXtract)](http://science.arm.gov/~cflynn/ARM_Tested_Tools/) is a scaled-down version of ANDX -- it is designed to only extract ASCII data. All features of ANDX pertaining to non-graphic data extraction are included in ANAX.

ANTS {#ANTS}
---------------------------

The ARM Program has developed [ANTS (ARM NetCDF Tool Suite)](http://science.arm.gov/~cflynn/ANTS/), a collection of netCDF tools and utilities providing various means of creating and modifying netcdf files. ANTS is based on nctools written by Chuck Denham. The utilities within nctools were modified to compile with version 3.5 of the netCDF library, the command syntax was modified for consistency with other tools, and changes were made to accommodate ARM standard netCDF.

The original functions from nctools were intended mainly for the creation, definition, and copying of fundamental netCDF elements. ARM added others which focus on manipulation of data within existing netCDF files. Additional functions have special support for multi-dimensional data such as "slicing" cross sections from multi-dimensional variable data or joining lesser-dimensional fields to form multi-dimensional structures. Functions have been added to support execution of arithmetic and logical operations, bundling or splitting netCDF files, comparing the structure or content of files, and so on.

Essentially every type of netCDF library function call is exercised in ANTS. In this way then, this open-source collection of tools also represents a library of coding examples for fundamental netCDF tasks. See the [website](http://science.arm.gov/~cflynn/ANTS/) for more information.

ARGOS {#ARGOS}
-----------------------------

[ARGOS](http://www.lapeth.ethz.ch/argos/index.html) (interActive thRee-dimensional Graphics ObServatory) is a new IDL-based interactive 3D visualization tool, developed by [David N. Bresch](http://www.lapeth.ethz.ch/~david/index.html) and [Mark A. Liniger](http://www.lapeth.ethz.ch/~mark/index.html) at the Institute for Atmospheric Science at the Swiss Federal Institute of Technology, ETH, ZÃ¼rich.

A highly optimized graphical user interface allows quick and elegant creation of even complex 3D graphics (volume rendering, isosurfaces,...), including Z-buffered overlays (with hidden lines), light and data shading, Xray images, 3D trajectories, animations and virtual flights around your data, all documented in a full on-line [html-help](http://www.lapeth.ethz.ch/argos/argos_general.html). The netCDF data format is preferred, but any other format can be read by providing an IDL (or FORTRAN or C or C++) interface. Some toolboxes (for atmospheric model output, trajectory display, radar data) have already been written, others might easily be added (in IDL, FORTRAN or C code). All interactive activities are tracked in a script, allowing quick reconstruction of anything done as well as running ARGOS in batch script mode.

Information about [copyright and licensing conditions](http://www.lapeth.ethz.ch/argos/argos_copyright.html) are available. For further information and installation, please E-mail to: bresch@atmos.umnw.ethz.ch

CDAT {#CDAT}
---------------------------

The [Climate Data Analysis Tool (CDAT)](http://cdat.sf.net), developed
by the [Program for Climate Model Diagnosis and Intercomparison
(PCMDI)](http://www-pcmdi.llnl.gov/) at Lawrence Livermore National
Laboratory, provides the capabilities needed to analyze model data,
perform complex mathematical calculations, and graphically display the
results. It provides the necessary tools to diagnose, validate, and
intercompare large observational and global climate model data sets.
It includes the ability to ingest large climate datasets in netCDF, HDF,
DRS, and GrADS/GRIB format; the Visualization and Computation System
(VCS) module, visually displays and animates ingested or created data;
and the Library of AMIP Data Transmission Standards (LATS) module
outputs data in the machine-independent netCDF or GrADS/GRIB file
formats.

In addition, the Command Line Interface (CLI) module allows CDAT to
receive argument and function input via the command line, and the
Graphical User Interface (GUI) allows CDAT to receive argument and
function input via a point-and-click environment.

The software, which runs as a standalone process or within PCMDI's
Visualization and Computation System (VCS), provides climate scientists
with an easy and fast method to read different file formats, and to
analyze and graphically display climate data in an integrated fashion.
CDAT includes a set of pre-defined functions to allow the user to
manipulate the data and send the output to a file which can be viewed as
an image, or as a collection of images in an animation. The software has
a gradual learning curve, allowing the novice user to quickly obtain
useful results.

CDFconvert {#CDFconvert}
---------------------------------------

The [MRG CDFconvert
package](http://www.atmos.albany.edu/facstaff/rmctc/cdf_cvt/) provided
by the Mesoscale Research Group, McGill University/SUNY Albany, is
designed to address data conversion issues for gridded datasets stored
under the
[COARDS](http://ferret.wrc.noaa.gov/noaa_coop/coop_cdf_profile.html)
convention. CDFconvert converts regular Cylindrical Equidistant
(Lat/Long) and Gaussian (Spherical) netCDF grids into either the
Canadian [RPN Standard
File](http://www.cmc.ec.gc.ca/rpn/modcom/si/libraries/rmnlib/fstd/index.html)
or [GEMPAK](/software/gempak/index.html) file formats. MRG CDFconvert
has the flexibility to handle netCDF files generated by a number of
sources, including NCEP and ECMWF. User-definable conversion tables make
the extension of the package to different datasets possible.

cdfsync {#cdfsync}
---------------------------------

Joe Sirott of NOAA's Pacific Marine Environmental Laboratory has
developed cdfsync, a program that allows users to rapidly synchronize a
set of netCDF files over a network. Fast synchronization times are
achieved by only transmitting the differences between files. It is built
on the Open Source [rsync](http://samba.anu.edu.au/rsync/) program, but
contains a number of optimizations including:

-   Special handling of netCDF files for faster synchronization
    calculations
-   Much faster updates of large numbers of small netCDF files
-   In-place updates of large netCDF files

The latest version should run on Linux variants and Solaris.

More information is available at the [cdfsync
website](http://www.epic.noaa.gov/epic/software/cdfsync/).

CDO (Climate Data Operators) {#CDO}
--------------------------------------------------

Uwe Schulzweida at the Max Planck Institute for Meteorology has
developed [CDO](http://code.zmaw.de/projects/cdo), a collection of
Operators to manipulate and analyze Climate Data files. Supported file
formats include netCDF and GRIB. There are more than 350 operators
available. The following table provides a brief overview of the main
categories.

-   File information (info, sinfo, diff, ...)
-   File operations (copy, cat, merge, split\*, ...)
-   Selection (selcode, selvar, sellevel, seltimestep, ...)
-   Missing values (setctomiss, setmisstoc, setrtomiss)
-   Arithmetic (add, sub, mul, div, ...)
-   Mathematical functions (sqrt, exp, log, sin, cos, ...)
-   Comparison (eq, ne, le, lt, ge, gt, ...)
-   Conditions (ifthen, ifnotthen, ifthenc, ifnotthenc)
-   Field statistics (fldsum, fldavg, fldstd, fldmin, fldmax, ...)
-   Vertical statistics (vertsum, vertavg, vertstd, vertmin, ...)
-   Time range statistics (timavg, yearavg, monavg, dayavg, ...)
-   Field interpolation (remapbil, remapcon, remapdis, ...)
-   Vertical interpolation (ml2pl, ml2hl)
-   Time interpolation (inttime, intyear)

As an example of use of CDO, converting from GRIB to netCDF can be as
simple as

        cdo -f nc copy file.grb file.nc

or with relative time axis (for usage with GrADS)
        cdo -r -f nc copy file.grb file.nc

or using ECMWF reanalysis on a reduced grid
        cdo -R -f nc copy file.grb file.nc

More information is available on the [CDO
homepage](http://code.zmaw.de/projects/cdo).

CIDS Tools {#CIDS_Tools}
---------------------------------------

The Center for Clouds Chemistry and Climate
([C4](http://www-c4.ucsd.edu/)) Integrated Data Systems
([CIDS](http://www-c4.ucsd.edu/~cids/)) group has developed several
useful netCDF utilities:
-   cdf2idl: Writes an IDL script to read a NetCDF file.
-   cdf2c: Writes C code to read a NetCDF file.
-   cdf2fortran: Writes FORTRAN source code to read a NetCDF file.
-   cdf2asc: Dumps NetCDF data to an ASCII file.

The source for these utilities can be downloaded from [CIDS NetCDF
Visualization Tools
site](http://www-c4.ucsd.edu/~cids/software/visual.html).

CSIRO MATLAB/netCDF interface {#CSIRO-MATLAB}
------------------------------------------------------------

The [CSIRO MATLAB/netCDF interface](http://www.marine.csiro.au/sw/matlab-netcdf.html) is now
available from the [CSIRO Marine Laboratories](http://www.marine.csiro.au).
The CSIRO MATLAB/netCDF interface is run from within MATLAB and has a
simple syntax. It has options for automatically handling missing values,
scale factors, and permutation of hyperslabs. It is, however, limited to
retrieving data from, and information about, existing netCDF files.

The basis of the interface is a machine-dependent mex-file called
mexcdf53. Rather than call the mex-file directly users are advised to
employ both [Chuck Denham's netCDF toolbox](#NC4ML5) and the CSIRO
MATLAB/netCDF interface described here. For read-only access to existing
netCDF data, the CSIRO interface has a simpler syntax than the netCDF
Toolbox, but the latter may also be used to create and manipulate netCDF
variables and datasets.

EPIC {#EPIC}
---------------------------

NOAA's Pacific Marine Environmental Laboratory
([PMEL](http://www.pmel.noaa.gov/)) has developed the
[EPIC](http://www.pmel.noaa.gov/epic/) software package for
oceanographic data. EPIC provides graphical display and data field
manipulation for multi-dimensional netCDF files (up to 4 dimensions).
PMEL has been using this software on Unix and VMS several years. At
present, they have:

-   a data file I/O library (
    [epslib](http://www.pmel.noaa.gov/epic/eps-manual/epslib_toc.html),
    which is layered on top of the netCDF library).
-   epslib allows transparent access to multiple data file formats
-   a [MATLAB MexEPS
    interface](http://www.epic.noaa.gov/epic/software/mexeps.htm) for
    using any supported EPIC file with MATLAB
-   [suite of EPIC
    programs](http://www.epic.noaa.gov/epic/software/ep_programs.htm)
    for graphics and analysis of hydrographic profile data and time
    series data.

This software was developed on Sun/Unix and is also supported for
DEC/Ultrix and VAX/VMS as a system for data management, display and
analysis system for observational oceanographic time series and
hydrographic data. The EPIC software includes over 50 programs for
oceanographic display and analysis, as well as utilities for putting
in-situ or observational data on-line (with on-the-fly graphics and data
download) on the WWW.
The developers are interested in coordinating with others who may be
developing oceanographic software for use with netCDF files. The EPIC
software is available via anonymous FTP from ftp.noaapmel.gov in the
epic/ and /eps directories. To obtain the EPIC software, please see Web
pages at <http://www.pmel.noaa.gov/epic/download/index.html>. For
information about EPIC, please see the Web pages at
<http://www.pmel.noaa.gov/epic/index.html>. Contact epic@pmel.noaa.gov,
or Nancy Soreide, nns@noaapmel.gov, for more information.

Excel Use
------------------------------------

Several packages are available for accessing netCDF data from Microsoft
Excel, including the [netcdf4excel](#netcdf4excel) add-in for Excel, and
a [Scientific Dataset (SDS) Library](#SDS) that supports a DataSetEditor
add-in for Excel to view and modify various forms of data, including
netCDF.

EzGet {#EzGet}
-----------------------------

A FORTRAN library called
[EzGet](http://www-pcmdi.llnl.gov/ktaylor/ezget/ezget.html) has been
developed at [PCMDI](http://www-pcmdi.llnl.gov/PCMDI.html) to facilitate
retrieval of modeled and observed climate data stored in popular formats
including [DRS](http://www-pcmdi.llnl.gov/drach/DRS.html),
[netCDF](/software/netcdf/), [GrADS](http://grads.iges.org/grads), and,
if a control file is supplied,
[GRIB](ftp://nic.fb4.noaa.gov/pub/nws/nmc/docs/gribed1/). You can
specify how the data should be structured and whether it should undergo
a grid transformation before you receive it, even when you know little
about the original structure of the stored data (e.g., its original
dimension order, grid, and domain).
The EzGet library comprises a set of subroutines that can be linked to
any FORTRAN program. EzGet reads files through the
[cdunif](http://www-pcmdi.llnl.gov/drach/cdunif.html) interface, but use
of EzGet does not require familiarity with cdunif. The main advantages
of using EzGet instead of the lower level cdunif library include:

-   Substantial error trapping capabilities and detailed error messages
-   Versatile capability of conveniently selecting data from specified
    regions (e.g., oceans, North America, all land areas north of 45
    degrees latitude, etc.)
-   Ability to map data to a new grid at the time it is retrieved by
    EzGet
-   Automatic creation of \`\`weights'' for use in subsequent averaging
    or masking of data
-   Increased control in specifying the domain of the data to be
    retrieved.

For more information about EzGet, including instructions for downloading
the documentation or software, see the EzGet home page at
<http://www-pcmdi.llnl.gov/ktaylor/ezget/ezget.html>. For questions or
comments on EzGet, contact Karl Taylor (taylor13@llnl.gov).

FAN
-------------------------

[FAN (File Array Notation)](fan_utils.html) is Harvey
Davies' package for extracting and manipulating array data from netCDF
files. The package includes the three utilities nc2text, text2nc, and
ncrob for printing selected data from netCDF arrays, copying ASCII data
into netCDF arrays, and performing various operations (sum, mean, max,
min, product, ...) on netCDF arrays. A library (fanlib) is also included
that supports the use of FAN from C programs. The package is available
via anonymous FTP from
<ftp://ftp.unidata.ucar.edu/pub/netcdf/contrib/fan.tar.Z>. Questions and
comments may be sent to Harvey Davies, harvey.davies@csiro.au.

FERRET {#FERRET}
-------------------------------

[FERRET](http://ferret.wrc.noaa.gov/Ferret/) is an interactive computer
visualization and analysis environment designed to meet the needs of
oceanographers and meteorologists analyzing large and complex gridded
data sets. It is available by anonymous ftp from abyss.pmel.noaa.gov for
a number of computer systems: SUN (Solaris and SUNOS), DECstation
(Ultrix and OSF/1), SGI, VAX/VMS and Macintosh (limited support), and
IBM RS-6000 (soon to be released).
FERRET offers a Mathematica-like approach to analysis; new variables may
be defined interactively as mathematical expressions involving data set
variables. Calculations may be applied over arbitrarily shaped regions.
Fully documented graphics are produced with a single command. Graphics
styles included line plots, scatter plots, contour plots, color-filled
contour plots, vector plots, wire frame plots, etc. Detailed controls
over plot characteristics, page layout and overlays are provided. NetCDF
is supported both as an input and an output format.

Many excellent software packages have been developed recently for
scientific visualization. The features that make FERRET distinctive
among these packages are Mathematica-like flexibility, geophysical
formatting (latitude/longitude/date), "intelligent" connection to its
data base, special memory management for very large calculations, and
symmetrical processing in 4 dimensions. Contact Steve Hankin,
hankin@noaapmel.gov, for more information.

Fimex {#fimex}
-----------------------------

Heiko Klein (Norwegian Meteorological Institute) has developed the
[fimex](https://wiki.met.no/fimex/start) (File Interpolation,
Manipulation, and EXtraction) C++ library for gridded geospatial data.
It converts between several data formats (currently netCDF, NcML, GRIB1
or GRIB2, and felt). Fimex also enables you to change the projection and
interpolation of scalar and vector grids, to subset the gridded data,
and to extract only parts of the files. Fimex supports a growing list of
other [features](https://wiki.met.no/fimex/features), including support
for most NcML features and for netCDF-4 compression.

For simple usage, Fimex also comes with the command line program fimex.

Documentation and downloads are available from the [fimex web
site](http://wiki.met.no/fimex/).

FWTools (GIS Binary Kit for Windows and Linux) {#fwtools}
------------------------------------------------------------------------

[FWTools](http://fwtools.maptools.org/) is Frank Warmerdam's set of Open
Source GIS binaries for Windows (win32) and Linux (x86 32bit) systems.
The kits are intended to be easy for end users to install and get going
with, and include OpenEV, GDAL, MapServer, PROJ.4 and OGDI as well as
some supporting components. FWTools aims to track the latest development
versions of the packages included as opposed to official releases, "to
give folks a chance to use the *latest and greatest*".

GDAL {#GDAL}
---------------------------

Frank Warmerdam's [GDAL](http://www.remotesensing.org/gdal/index.html)
is a translator library for raster geospatial data formats that is
released under an X/MIT style Open Source license. As a library, it
presents a [single abstract data
model](http://www.remotesensing.org/gdal/gdal_datamodel.html) to the
calling application for all supported formats. The related
[OGR](http://www.remotesensing.org/gdal/ogr) library (which lives within
the GDAL source tree) provides a similar capability for simple features
vector data.

GDAL is in active use in several projects, and includes roughly 40
format drivers, including a translator for netCDF (read/write). Other
translators include GeoTIFF (read/write), Erdas Imagine (read/write),
ESRI .BIL (read), .aux labeled raw (read/write), DTED (read), SDTS DEM
(read), CEOS (read), JPEG (read/write), PNG (read/write), Geosoft GXF
(read) and Arc/Info Binary Grid (read). A full list is available in
[Supported
Formats](http://www.remotesensing.org/gdal/formats_list.html).

GDAL has recently included support for the netCDF-4 enhanced data model
and netCDF-4 format, as well as improved support for recent additions to
the CF conventions.

As an example of the use of GDAL, converting an ArcInfo ASCII grid to
netCDF (GMT conventions) as easy as:

       gdal_translate arc_ascii.grd -of GMT gmt_grid.nc

GDL (GNU Data Language) {#GDL}
---------------------------------------------

[GDL](http://gnudatalanguage.sourceforge.net/) is a free implementation
of most of the programming language supported by [IDL](#IDL)
(Interactive Data Language). GDL supports the netCDF-3 API.

Gfdnavi (Geophysical fluid data navigator) {#Gfdnavi}
--------------------------------------------------------------------

[Gfdnavi](http://www.gfd-dennou.org/arch/davis/gfdnavi/index.en.htm) is
a web-based tool to archive, share, distribute, analyze, and visualize
geophysical fluid data and knowledge. The software is under development
by members of the GFD Dennou Club, including T. Horinouchi (RISH, Kyoto
U.), S. Nishizawa (RIMS, Kyoto U.), and colleagues. Gfdnavi uses a
metadata database for managing and analyzing data and visualizations. It
also permits publishing data for web access and will soon support access
to data on other Gfdnavi servers. Web service APIs are now under
development. A presentation [Introducing
Gfdnavi](http://www.gfd-dennou.org/arch/davis/gfdnavi/presen/2007-03-05_GfdnaviIntro.En/pub/)
describes the architecture and shows examples of use.

Gfdnavi is dependent on two technologies:

-   [Ruby on Rails](http://www.rubyonrails.com/), a framework for web
    applications, and
-   [the Dennou Ruby Project](http://ruby.gfd-dennou.org/), a collection
    of tools for geophysical data. These tools include
    [GPhys](http://ruby.gfd-dennou.org/products/gphys/) software to
    handle GRIB, GrADS, and netCDF data uniformly.

As an example of this technology, Takuji Kubota has established [a
Gfdnavi server](http://www.gsmap.aero.osakafu-u.ac.jp/gfdnavi/) for the
Global Satellite Mapping of Precipitation
([GSMaP](http://www.radar.aero.osakafu-u.ac.jp/~gsmap/index_english.html))
project.

GMT {#GMT}
-------------------------

[GMT](http://gmt.soest.hawaii.edu/) (Generic Mapping Tools) is an open
source collection of about 60 tools for manipulating geographic and
Cartesian data sets (including filtering, trend fitting, gridding,
projecting, etc.) and producing Encapsulated PostScript File (EPS)
illustrations ranging from simple x-y plots via contour maps to
artificially illuminated surfaces and 3-D perspective views. GMT
supports 30 map projections and transformations and comes with support
data such as coastlines, rivers, and political boundaries. GMT is
developed and maintained by Paul Wessel and Walter H. F. Smith with help
from a global set of volunteers, and is supported by the National
Science Foundation. It is released under the GNU General Public License.

The package can access COARDS-compliant netCDF grids as well as ASCII,
native binary, or user-defined formats. The GMT package is available via
anonymous ftp from several servers; see
[gmt.soest.hawaii.edu](http://gmt.soest.hawaii.edu) for installation
information.

Grace {#Grace}
-----------------------------

[Grace](http://plasma-gate.weizmann.ac.il/Grace/) is a tool to make
two-dimensional plots of scientific data, including 1D netCDF variables.
It runs under the X Window System and OSF Motif (recent versions of
LessTif are, by and large, fine, too). Grace runs on practically any
version of Unix. As well, it has been successfully ported to VMS, OS/2
and Win9\*/NT (some functionality may be missing, though). Grace is a
descendant of ACE/gr.
A few features of Grace are:

-   User defined scaling, tick marks, labels, symbols, line styles,
    colors.
-   Batch mode for unattended plotting.
-   Read and write parameters used during a session.
-   Regressions, splines, running averages, DFT/FFT,
    cross/auto-correlation, ...
-   Support for dynamic module loading.
-   Hardcopy support for PostScript, PDF, GIF, and PNM formats.
-   Device-independent Type1 font rastering.
-   Ability to read or write netCDF data.

GrADS {#GrADS}
-----------------------------

[GrADS](http://grads.iges.org/grads/grads.html) (Grid Analysis and
Display System) is an interactive desktop tool from
[COLA/IGES](http://grads.iges.org/cola.html) that is currently in use
worldwide for the analysis and display of earth science data. GrADS is
implemented on all commonly available UNIX workstations, Apple
Macintosh, and DOS or Linux based PCs, and is freely available via
anonymous ftp. GrADS provides an integrated environment for access,
manipulation, and display of earth science data in several forms,
including GRIB and netCDF. For more information, see the [GrADS User's
Guide](http://grads.iges.org/grads/gadoc/users.html).

Gri
-------------------------

Gri is an extensible plotting language for producing scientific graphs,
such as x-y plots, contour plots, and image plots. Dan Kelley of
Dalhousie University is the author of Gri, which can read data from
netCDF files as well as ASCII and native binary data. For more
information on Gri, see the URL <http://gri.sourceforge.net/>.

GXSM {#GXSM}
---------------------------

The GXSM is the Gnome X Scanning Microscopy project, it is a bit more
than just a piece of software (the GXSM itself), there is full hardware
support for DSP cards including open source DSP software and a growing
set of SPM related electronics. For more information, see
<http://gxsm.sourceforge.net/>.

HDF interface {#HDF_interface}
---------------------------------------------

The National Center for Supercomputing Applications (NCSA) has added the
netCDF interface to their [Hierarchical Data Format
(HDF)](http://hdf.ncsa.uiuc.edu/) software. HDF is an extensible data
format for self-describing files. A substantial set of applications and
utilities based on HDF is available; these support raster-image
manipulation and display and browsing through multidimensional
scientific data. An implementation is now available that provides the
netCDF interface to HDF. With this software, it is possible to use the
netCDF calling interface to place data into an HDF file. The netCDF
calling interface has not changed and netCDF files stored in XDR format
are readable, so existing programs and data will still be usable
(although programs will need to be relinked to the new library). There
is currently no support for the mixing of HDF and netCDF structures. For
example, a raster image can exist in the same file as a netCDF object,
but you have to use the Raster Image interface to read the image and the
netCDF interface to read the netCDF object. The other HDF interfaces are
currently being modified to allow multi-file access, closer integration
with the netCDF interface will probably be delayed until the end of that
project.
Eventually, it will be possible to integrate netCDF objects with the
rest of the HDF tool suite. Such an integration will then allow tools
written for netCDF and tools written for HDF to both interact
intelligently with the new data files.

HDF-EOS to netCDF converter {#HDF-EOS}
-----------------------------------------------------

The Goddard Earth Sciences Data and Information Services Center ([GES
DISC](http://disc.gsfc.nasa.gov)) has developed an on-the-fly HDF-EOS to
netCDF/CF converter for the following products, making them easier to
use in the [Unidata IDV](#IDV) and
[McIDAS-V](http://www.ssec.wisc.edu/mcidas/software/v/):

-   AIRS Level 2 (scene) profiles of moisture, air temperature and trace
    gases
-   AIRS Level 3 (global grid) profiles of moisture, air temperature and
    trace gases
-   OMI UV-B at the surface
-   TOMS ozone and aerosols

[Instructions](http://disc.gsfc.nasa.gov/services/NetCDFConversionforIDVandMcIDAS-V.shtml)
are available for searching and converting these data. More information
on AIRS products is available at
<http://disc.gsfc.nasa.gov/AIRS/index.html>.

HIPHOP {#HIPHOP}
-------------------------------

[HIPHOP](http://www.knmi.nl/onderzk/atmosam/English/Service/hiphop/hiphop.html),
developed by Dominik Brunner, is a widget based IDL application that
largely facilitates the visualization and analysis of 2D, 3D, and 4D
atmospheric science data, in particular atmospheric tracer distributions
and meteorological fields.
Graphical output of (atmospheric model) data can be quickly generated in
a large number of different ways, including horizontal maps at selected
model or pressure levels, vertical north-south, east-west, or slant
cross-sections (including zonal averages), time slices, animations, etc.
It also allows mathematical operations on the existing fields to
generate new fields for further analysis, and it can be run as a batch
application.

The program handles data in netCDF, HDF and GRIB format. Interfaces to
other data formats (e.g. ASCII and binary data) can be added easily.

Beginning with Version 4.0, it also supports the ability to overlay
meteorological fields on a number of different satellite images, and to
draw air parcel trajectories.

Hyperslab OPerator Suite (HOPS) {#HOPS}
---------------------------------------------------------------------------------

Hyperslab OPerator Suite
([HOPS](http://www.cgd.ucar.edu/gds/svn/hyperslab.html)), developed by
R. Saravanan at NCAR, is a bilingual, multi-platform software package
for processing data in netCDF files conforming to the NCAR-CCM format or
the NCAR Ocean Model format. HOPS is implemented in [IDL](#IDL), the
widely-used commercial interpreted language, and also in
[Yorick](ftp://ftp-icf.llnl.gov/pub/Yorick/), a public-domain
interpreted language that is freely available from the Lawrence
Livermore National Laboratory. The IDL version of HOPS should run on any
platform supported by IDL. The Yorick version too runs on most common
UNIX platforms, such as Sun, SGI, Cray, and LINUX computers.
HOPS is not a monolithic program, but a suite of operators that act on
data units called "hyperslabs". The design of HOPS is object-oriented,
rather than procedure-oriented; the operators treat the numeric data and
the associated meta-data (like coordinate information) as a single
object.

Note that HOPS is not a general purpose netCDF utility and works only
for the NCAR CSM netCDF formats. For more information, check the [HOPS
home page](http://www.cgd.ucar.edu/gds/svn/hyperslab.html).

iCDF (imports chromatographic netCDF data into MATLAB) {#iCDF}
-----------------------------------------------------------------------------

Klavs M. Sørensen, Thomas Skov and Rasmus Bro (Faculty of Life Sciences,
University of Copenhagen) have developed
[iCDF](http://www.models.life.ku.dk/source/iCDF/index.asp), a free and
documented toolbox for importing chromatographic data in the
netCDF-based format that most manufacturers of chromatographic software
support.

The iCDF software is currently for XC-MS data (X: GC, LC, HPLC), but
soon it will be able to import data using other detectors as well. It
can be used to open netCDF files from many different instruments (e.g.
Agilent, Bruker) and many chromatographic software packages (e.g.
ChemStation).

For more information, see the paper

> Skov T and Bro R. (2008) Solving fundamental problems in
> chromatographic analysis Analytical and Bioanalytical Chemistry, 390
> (1): 281-285.

IDV (Integrated Data Viewer) {#IDV}
--------------------------------------------------

Unidata's [Integrated Data Viewer (IDV)](/software/idv/) is a Java
application (for Java 1.4 or later) that can be used to display a
variety of netCDF files, particularly well formatted, geolocated
datasets. Features include:

-   Access to local and remote netCDF files and a variety of [other data
    formats](/software/idv/docs/userguide/data/DataSources.html)
-   Slicing and probing of multidimensional data
-   Support for netCDF conventions (CF, COARDS, NUWG, AWIPS)
-   InstallAnywhere installers for easy download and installation
-   Save display state to a bundle for easy recreation of views
-   Support for non-gridded data through the [Common Data Model
    (CDM)](/software/netcdf-java/CDM/)

The IDV uses the [VisAD Java
library](http://www.ssec.wisc.edu/~billh/visad.html) for interactive and
collaborative visualization and analysis and the [netCDF Java
library](/software/netcdf-java/) for reading and manipulating netCDF
files.

Ingrid {#Ingrid}
-------------------------------

[Ingrid](http://ingrid.ldgo.columbia.edu/), by M. Benno Blumenthal
\<benno@ldeo.columbia.edu\>, is designed to manipulate large datasets
and model input/output. It can read data from its data catalog, a netCDF
file, or a directly attached model, and output the data, either by
feeding it to a model, creating a netCDF file, or creating plots and
other representations of the data.

Ingrid has a number of filters which allow simple data manipulations,
such as adding two datasets together, smoothing, averaging, and
regridding to a new coordinate. In addition to netCDF, it also reads
HDF, CDF, VOGL, and SGI GL.

Ingrid is currently running as a WWW daemon that can be accessed through
<http://rainbow.ldgo.columbia.edu/datacatalog.html> to see some of its
capabilities on a climate data catalog maintained by the [Climate
Group](http://rainbow.ldeo.columbia.edu/) of the [Lamont-Doherty Earth
Observatory](http://www.ldeo.columbia.edu/) of Columbia University. To
quote the introduction:

> The Data Catalog is both a catalog and a library of datasets, i.e. it
> both helps you figure out which data you want, and helps you work with
> the data. The interface allows you to make plots, tables, and files
> from any dataset, its subsets, or processed versions thereof.
>
> This data server is designed to make data accessible to people using
> WWW clients (viewers) and to serve as a data resource for WWW
> documents. Since most documents cannot use raw data, the server is
> able to deliver the data in a variety of ways: as data files (netCDF
> and HDF), as tables (html), and in a variety of plots (line, contour,
> color, vector) and plot formats (PostScript and gif). Processing of
> the data, particularly averaging, can be requested as well.
>
> The Data Viewer in particular demonstrates the power of the Ingrid
> daemon.

Ingrid currently runs on Linux, for which binaries are available. CVS
access to the current source can be arranged.

 Intel Array Visualizer {#IntelArrayVisualizer}
--------------------------------------------------------------

The [Intel® Array
Visualizer](http://www.intel.com/cd/software/products/asmo-na/eng/compilers/226277.htm)
and Intel® Array Viewer are available as [free
downloads](http://www.intel.com/cd/software/products/asmo-na/eng/compilers/226277.htm)
for Windows platforms. They offer an application and a set of software
tools and components, which include C, Fortran, and .Net libraries, for
developing scientific visualization applications and for creating
interactive graphs of array data in various formats, including HDF and
netCDF.

IVE {#IVE}
-------------------------

[IVE (Interactive Visualization
Environment)](http://www.atmos.washington.edu/ive/) is a software
package designed to interactively display and analyze gridded data. IVE
assumes the data to be displayed are contained in one- two-, three- or
four-dimensional arrays. By default, the numbers within these arrays are
assumed to represent grid point values of some field variable (such as
pressure) on a rectangular evenly spaced grid. IVE is, nevertheless,
capable of displaying data on arbitrary curvilinear grids.

If the data points are not evenly spaced on a rectangular grid, IVE must
be informed of the grid structure, either by specifying "attributes" in
the data input or by specifying the coordinate transform in a user
supplied subroutine. Stretched rectangular grids (which occur when the
stretching along a given coordinate is a function only of the value of
that coordinate) can be accommodated by specifying one-dimensional
arrays containing the grid-point locations along the stretched
coordinate as part of the IVE input data. Staggered meshes can also be
accommodated by setting "attributes" in the input data. The structure of
more complicated curvilinear grids must be communicated to IVE via user
supplied "transforms," which define the mapping between physical space
and the array indices.

Since four-dimensional data cannot be directly displayed on a flat
computer screen, it is necessary to reduced the dimensionality of the
data before it is displayed. One of IVE's primary capabilities involves
dimension reduction or "data slicing." IVE allows the user to display
lower-dimensional subsets of the data by fixing a coordinate or by
averaging over the coordinate.

IVE currently has the capability to display

-   scalar fields as
    -   2D scalar plots
    -   1D scalar plots
    -   vertical soundings
    -   a single point value
-   vector fields as 2D vector plots

IVE lets you overlay plots, loop plots, and control a wide variety of
display parameters.

IVE also can perform algebraic computations on the gridded data and can
calculate derivatives. More complicated computations can be performed in
user supplied subroutines.

IVE uses NetCDF for the data input format, and uses the [NCAR Graphics
Library](http://ngwww.ucar.edu/ng/) to produce graphical output. IVE is
[available](http://www.atmos.washington.edu/ive/getting.html) as source
via anonymous ftp; and as binary on request for licensees of NCAR
graphics.

JSON format with the ncdump-json utility {#JSON}
---------------------------------------------------------------

Josep Llodrà has developed a program to output the contents of a
netCDF-3 or netCDF-4 file in JSON (JavaScript Object Notation). It is
based on Unidata's NCDUMP utility, and it keeps the original ncdump
functionality, unless the "-j" option is used to specify JSON output.

The program and source are available from
<https://github.com/jllodra/ncdump-json> .

Java interface {#Java_interface}
-----------------------------------------------

The [NetCDF-Java 4.2 Library](/software/netcdf-java/) is a Java
interface to netCDF files, as well as to many other types of scientific
data formats. It is freely available and the source code is released
under the (MIT-style) netCDF C library license. Previous versions use
the GNU Lesser General Public License (LGPL).

The library implements a Common Data Model
([CDM](/software/netcdf-java/CDM/)), a generalization of the netCDF,
OpenDAP and HDF5 data models. The library is a prototype for the
netCDF-4 project, which provides a C language API for the "data access
layer" of the CDM, on top of the HDF5 file format. The NetCDF-Java
library is a 100% Java framework for *reading* netCDF and other file
formats into the CDM, as well as *writing* to the netCDF-3 file format.
The library also implements
[NcML](http://www.unidata.ucar.edu/software/netcdf/ncml/), which allows
you to add metadata to CDM datasets, as well as to create virtual
datasets through aggregation.

Kst (2D plotting tool) {#KST}
--------------------------------------------

[Kst](http://kst-plot.kde.org) is an open-source, cross-platform 2D
plotting tool focused on performance and ease of use. Packages for
Windows, various Linux distributions and Mac OS X are
[available](http://sourceforge.net/projects/kst/files/), as well as the
complete source code and CMake-based build files. A more detailed
presentation of Kst can be found on the web page at
<http://kst-plot.kde.org>, including numerous screenshots and all the
useful download links.

Kst is characterized by the following features:

-   Outstanding performance: curves with millions of points are no
    problem
-   Plotting of live streams
-   Out-of-the box support for a variety of formats (currently ASCII,
    netCDF, dirfile, Qimage-supported types, fits images)
-   User-friendly with a modern and consistent user interface
-   A set of unique tools to boost efficiency, including a data import
    wizard, capacity to edit multiple objects at once or the "Change
    Data File" tool to compare multiple experiments easily
-   An active community
-   Easily expandable for new data formats or data analysis algorithms
    thanks to a plugin-based architecture
-   Available on Windows, Linux, and Mac OSX

Labview interface {#Labview-API}
-----------------------------------------------

A netCDF Labview interface, implemented in the Labview programming
language is available. The software includes A graphical user interface
for editing netCDF data and conversion to other data formats. The
package was developed and is maintained by L. F. Hwang of Sun Yat-sen
University in China. For more information and to download the source
code, see the [NetCDFLabview web
site](https://sourceforge.net/projects/netcdflabview/).

MBDyn (MultiBody Dynamics) #{MBDyn}
--------------------------------------------------

[MBDyn](http://www.aero.polimi.it/~mbdyn/) is an open-source MultiBody
Dynamics analysis system developed at the Dipartimento di Ingegneria
Aerospaziale of the University "Politecnico di Milano", Italy. It uses
netCDF as its primary output format.

MBDyn features the integrated multidisciplinary analysis of multibody,
multiphysics systems, including nonlinear mechanics of rigid and
flexible constrained bodies, smart materials, electric networks, active
control, hydraulic networks, essential fixed-wing and rotorcraft
aerodynamics. It allows users to simulate the behavior of heterogeneous
mechanical, aero-servo-elastic systems based on first principles
equations. It is being actively developed and used in the aerospace and
automotive fields for dynamics analysis and simulation of complex
systems. Dynamic linking of user-defined modules is heavily exploited to
let users extend the feature library.

Max_diff_nc ${Maxdiffnc}
-------------------------------------------

This is a program which compares two NetCDF files. Variables with the
same ID in the two files are assumed to be of the same type and have the
same shape. For each such couple of variables, the program computes the
maximum of the absolute value of the difference, and the maximum of the
absolute value of the relative difference. The program also tells you at
what location (the subscript list of the array) the maximum difference
is reached.

The web page for this program is:
<http://web.lmd.jussieu.fr/~lglmd/Max_diff_nc>

This is a freely available tool.

MeteoExplorer {#MeteoExplorer}
---------------------------------------------

[MeteoExplorer](http://www.eastmodelsoft.com/index_en.htm), developed by
Lianqing Yu at China Meteorological Administration, is a cross-platform
software application for analyzing and rendering atmospheric science and
geoscience data. It supports popular data formats including WMO
GRIB1/GRIB2, NetCDF, and MICAPS, and provides basic GIS functionalities.
Developed with C++, Meteo Explorer targets multiple computing platforms
including Microsoft Windows, GNU Linux, and SGI IRIX operating systems.

The primary features include:

-   Graphics layer management (navigation and animation)
-   Objective analysis of physical elements in surface or upperair
    soundings data
-   Isoline analysis and shading of grid field
-   Streamline analysis of wind field
-   Computation of physics elements
-   NetCDF data process and display
-   GRIB1/GRIB2 data process and display
-   MICAPS data process and display
-   Satellite nephogram data display and animation, support AWX, GPF and
    HDF format
-   Interactive composition of synoptic chart (command undo/redo,
    automatic save)
-   Map zoom, pan, projection and clipping
-   Full screen display and zoom to area
-   Quick navigation via thumbnail view of graphics layers
-   Save screen shot as image file (support formats: BMP, JPG, PNG)
-   Vector graphics exported to clipboard or saved as EMF file (Windows
    version only)
-   Remote desktop connection support
-   System configuration (dynamic menu)
-   Fast switch of user interface language on the fly

For more information, please visit [MeteoExplorer's home
page](http://www.eastmodelsoft.com/software/mexplorer.htm) or contact
the support staff via meteoexplorer@hotmail.com .

MeteoInfo {#MeteoInfo}
-------------------------------------

For better cross-platform support,
[MeteoInfo](http://www.meteothinker.com) has recently been re-developed
using Unidata's NetCDF Java library. MeteoInfo is GIS software for
visualization and analysis of spatial and meteorological data. The Java
edition can be run in Windows, Mac OS, Linux, and Unix systems. The
Groovy script engine was coupled in the software, so users can write
Groovy script to run the software automatically for analysis with
complex steps.

Download: <http://www.meteothinker.com/>

Java 6 is needed to run the software.

MexEPS {#MexEPS}
-------------------------------

[PMEL](http://www.pmel.noaa.gov/) has developed a MATLAB interface,
[MexEPS](http://www.epic.noaa.gov/epic/software/mexeps.htm), which
supports several netCDF file conventions, including [those adopted by
PMEL](ftp://ftp.unidata.ucar.edu/pub/netcdf/Conventions/PMEL-EPIC/).
Many styles of time axes are supported and time manipulation routines
ease the use of the time axis in MATLAB. The MexEPS package supports the
following data formats:
-   reading, writing and editing netCDF files;
-   reading and writing Classic EPIC files
-   reading formatted ASCII files

It includes:
-   VARIABLE, AXIS, ATTRIBUTE manipulation routines
-   TIME manipulation
    -   TIME enters MATLAB as YYMMDDhhmmss.fff
    -   Can be converted to netCDF udunits time convention (e.g. days
        *since* 1990-01-01 00:00:00)
-   [MATLAB help](ftp://ftp.pmel.noaa.gov/eps/mexeps/help-m/) and
    [example scripts](ftp://ftp.pmel.noaa.gov/eps/mexeps/examples/)
    using MexEPS
-   **ASCII2MAT** mexFunction, which reads a formatted file into MATLAB
    as a matrix

The MexEPS package is freely available in PMEL's anonymous ftp directory
<ftp://ftp.pmel.noaa.gov/eps/mexeps/>

If you have any questions or comments, please contact the author, Willa
Zhu [(willa@pmel.noaa.gov)](mailto:willa@pmel.noaa.gov) or Nancy Soreide
(nns@pmel.noaa.gov).

MEXNC and SNCTOOLS {#MEXNC}
------------------------------------------

John Evans of Rutgers University maintains MEXNC and developed SNCTOOLS.
[MEXNC](http://mexcdf.sourceforge.net/) is a mexfile interface to NetCDF
files for MATLAB that has roughly a one-to-one equivalence with the C
API for netCDF.
[SNCTOOLS](http://mexcdf.sourceforge.net/tutorial/index.html) is a set
of higher-level m-files that sit atop MEXNC, shielding the user from
such low level netCDF details as file IDs, variable IDs, and dimension
IDs. The general philosophy behind SNCTOOLS is providing the ability to
read and write data without trying to invent a new syntax.

Mirone (Windows MATLAB-based display) {#Mirone}
--------------------------------------------------------------

Joaquim Luis of Universidade do Algarve has developed
[Mirone](http://w3.ualg.pt/~jluis/mirone/), a Windows MATLAB-based
framework tool that allows the display and manipulation of a large
number of grid/images formats through its interface with the
[GDAL](http://remotesensing.org/gdal/) library. Its main purpose is to
provide users with an easy-to-use graphical interface to manipulate
[GMT](http://gmt.soest.hawaii.edu/) grids. In addition it offers a wide
range of tools dedicated to topics in the earth sciences, including
tools for multibeam mission planning, elastic deformation studies,
tsunami propagation modeling, earth magnetic field computations and
magnetic Parker inversions, Euler rotations and poles computations,
plate tectonic reconstructions, and seismicity and focal mechanism
plotting. The high quality mapping and cartographic capabilities for
which GMT is renowned is guaranteed through Mirone's ability to
automatically generate GMT cshell scripts and dos batch files.

Although Mirone is written in MATLAB, a stand-alone version to run under
Windows is also provided. Regrettably this version is not as efficient
as the native MATLAB code but provides a solution for users that don't
have MATLAB.

Also see\
 J. F. Luis. Mirone: A multi-purpose tool for exploring grid data.
Computers & Geosciences, 33, 31-41, 2007.

ncBrowse {#ncBrowse}
-----------------------------------

Donald Denbo of NOAA's Pacific Marine Environmental Laboratory has
developed and made available
[ncBrowse](http://www.epic.noaa.gov/java/ncBrowse), a Java application
(JDK1.2) that provides flexible, interactive graphical displays of data
and attributes from a wide range of netCDF data file conventions.
Features include:

-   Designed to work with arbitrary netCDF files.
-   Browses file using the EPIC and COARDS conventions.
-   Provides a "tree" view of the netCDF file.
-   Handles character variables.
-   Handles dimensions without an associated variable.
-   Uses sgt graphics to perform 1 and 2 dimensional cuts through data.
-   Save to file single variable as a "cdl" text file.
-   InstallAnywhere scripts for UNIX, Win32, and MacOS.
-   Currently uses Java 2 and Swing.

ncBrowse will run on any UNIX or Windows machine with a Java 2 (JDK1.2)
virtual machine installed. Automated installation scripts are available
for Windows and UNIX. Additional information on ncBrowse and download
instructions are available at <http://www.epic.noaa.gov/java/ncBrowse>.

Questions and suggestions should be directed to
\<[dwd@pmel.noaa.gov\>](mailto:dwd@pmel.noaa.gov). If you have problems
reading a netCDF file with ncBrowse, please send him a copy of the file
and he'll get ncBrowse to read it!

nccmp {#nccmp}
-----------------------------

Remik Ziemlinski of the NOAA Geophysical Fluid Dynamics Laboratory has
developed [nccmp](http://nccmp.sourceforge.net/), a tool to compare two
netCDF files. It can use MPI, include/exclude specific variables or
metadata and operates quickly. Highly recommended for regression testing
with large datasets. See the Web site <http://nccmp.sourceforge.net/>
for more information.

NCL {#NCL}
-------------------------

The [NCAR Command Language (NCL)](http://www.ncl.ucar.edu/) is an
intepreted programming language for scientific data analysis and
visualization developed and maintained in NCAR's [Computational and
Information Systems Laboratory](http://www.cisl.ucar.edu/).

NCL has many features common to modern programming languages, including
types, variables, operators, expressions, conditional statements, loops,
and functions and procedures. NCL also has features that are not found
in other programming languages, including those that handle the
manipulation of metadata, the configuration of visualizations, the
import of data from a variety of data formats, and an algebra that
supports array operations.

NCL has robust file input and output capabilities. It allows different
datasets of different formats (netCDF, netCDF-4 classic, HDF4, HDF4-EOS,
GRIB-1, and GRIB-2) to be imported into one uniform and consistent data
manipulation environment, which internally is the netCDF data format.
NCL doesn't place any restrictions or conventions on the organization of
input netCDF files.

NCL comes with many useful built-in functions and procedures for
processing and manipulating data. There are over 600 functions and
procedures that include routines for use specifically with climate and
model data, empirical orthogonal functions, Fourier coefficients,
wavelets, singular value decomposition, 1-, 2-, and 3-dimensional
interpolation, approximation, and regridding, and computer analysis of
scalar and vector global geophysical quantities.

The visualizations are publication-quality and highly customizable, with
hundreds of options available for tweaking the looks of your graphics.
NCL can generate contours, XY plots, vectors, streamlines, and can
overlay these plots on many different map projections. There are also
specialized functions for generating histograms, wind roses, meteograms,
skew-T plots, weather maps.

Included with the software are two command line tools: "ncl\_convert2nc"
for converting GRIB-1/2 or HDF files to netCDF files, and
"ncl\_filedump" which will dump the contents of a file format that NCL
recognizes (netCDF, GRIB-1/2, HDF, etc).

NCL is available under an open source license or in binary form for
several popular UNIX platforms, including (but not limited to) Linux,
MacOSX, and Windows/Cygwin.

Documentation and additional information on NCL are available from the
[NCL website](http://www.ncl.ucar.edu/), which contains hundreds of
[application examples](http://www.ncl.ucar.edu/Applications/) for one to
download. You can also contact Mary Haley, at <haley@ucar.edu> for more
information.

NCO {#NCO}
-------------------------

[NCO](http://nco.sourceforge.net) (netCDF operators) is a package of
command line operators that work on generic netCDF or HDF4 files:
-   ncap2 - arithmetic processor
-   ncatted - attribute editor
-   ncbo - binary operator
-   ncdiff - differencer
-   ncea - ensemble averager
-   ncecat - ensemble concatenator
-   ncflint - file interpolator
-   ncks - kitchen sink (extract, cut, paste, print data)
-   ncpdq - permute dimensions quickly
-   ncra - running averager
-   ncrcat - record concatenator
-   ncrename - renamer
-   ncwa - weighted averager

All operators may now be [OPeNDAP](http://www.opendap.org) clients. OPeNDAP
enables network transparent data access to any OPeNDAP server. Thus
OPeNDAP-enabled NCO can operate on remote files accessible through any
OPeNDAP server without transferring the files. Only the required data
(e.g., the variable or hyperslab specified) are transferred.

The source code is freely available from the [NCO home
page](http://nco.sourceforge.net/), as is the NCO User's Guide.

For more information, contact the author, Charlie Zender.

ncregrid {#ncregrid}
-----------------------------------

Patrick Jöckel of the Max Planck Institute for Chemistry has developed
**ncregrid**, a tool (written in FORTRAN-90) for data transfer of
gridded 2- and 3-dimensional (spatial) geophysical/geochemical scalar
fields between grids of different resolutions. The algorithm handles
data on rectangular latitude/longitude grids (not necessarily evenly
spaced) and vertical pressure hybrid grids of arbitrary resolution. The
input/output data format is netCDF. ncregrid is freely available without
any warranty under the GNU public license (GPL). ncregrid can be used as
a "stand-alone" program, and/or linked as an interface to a model, in
order to re-grid automatically the input from an arbitrary grid space
onto the required grid resolution.

More information is available on the web-page:
<http://www.mpch-mainz.mpg.de/~joeckel/ncregrid/index.html>.

nctoolbox (a MATLAB common data model interface) {#nctoolbox}
----------------------------------------------------------------------------

[nctoolbox](http://nctoolbox.github.io/nctoolbox/) is a MATLAB interface
that provides read-only access to [Common Data
Model](/software/netcdf-java/CDM/index.html) datasets. Under the hood,
nctoolbox uses Unidata's NetCDF-Java as the data access layer. This
allows nctoolbox to access to netCDF, OPeNDAP, HDF5, GRIB, GRIB2, HDF4,
and many (15+) other file formats and services using the same API. It
works with MATLAB 2008a and later. The nctoolbox software was developed
by Brian Schlining (MBARI), Rich Signell (USGS), Sachin Kumar Bhate
(freelance), and Alex Crosby (RPS/ASA).

ncdx {#ncdx}
---------------------------

Patrick Jöckel of the Max Planck Institute for Chemistry has developed
**ncdx**, a tool (written in FORTRAN-90) that scans a netCDF file and
makes it [OpenDX](#OpenDX) compliant. ncdx is freely available without
any warranty under the GNU public license (GPL). More information is
available on the web-page:
<http://www.mpch-mainz.mpg.de/~joeckel/ncdx/index.html>.

ncensemble {#ncensemble}
---------------------------------------

Alan Iwi, of Rutherford Appleton Laboratory, offers this command line
ensemble statistics utility. More information is available on the
web-page: <http://home.badc.rl.ac.uk/iwi/ncensemble/>.

ncview {#ncview}
-------------------------------

[Ncview](http://meteora.ucsd.edu/~pierce/ncview_home_page.html) is a
visual browser for netCDF files. Typically you would use ncview to get a
quick and easy, push-button look at your netCDF files. You can view
simple movies of the data, view along various dimensions, take a look at
the actual data values, change color maps, invert the data, etc. It runs
on UNIX platforms under X11, R4 or higher. For more information, check
out the [README](http://meteora.ucsd.edu/~pierce/docs/ncview.README)
file; you can also see a representative [screen
image](http://meteora.ucsd.edu/~pierce/docs/ncview.gif) (GIF, 66K) of
ncview in action.
The source may be downloaded from <ftp://cirrus.ucsd.edu/pub/ncview/>.
For more information, please contact the author, David W. Pierce at
<dpierce@ucsd.edu>.

netcdf4-js {#netcdf4-js}
-------------------------------
[netcdf4-js](https://www.npmjs.com/package/netcdf4) is a NodeJS addon for reading and writing the files in the Network Common Data Form (NetCDF) version <= 4, built upon the C-library for netcdf. It isavailable from npmjs at the link above, or directly from the [GitHub Repository](https://github.com/swillner/netcdf4-js).

NetCDF Toolbox for MATLAB-5 {#matlab5}
----------------------------------------------------

The [NetCDF Toolbox for MATLAB-5](http://mexcdf.sourceforge.net/),
originally developed by Charles R. Denham, combined netCDF-3 with
[MATLAB](http://www.mathworks.com/products/matlab/) to form an interface
that used MATLAB operator-syntax for arithmetic, logical, and
subscripting operations on netCDF entities. The NetCDF Toolbox is in
bug-fix-only mode, and is maintained by John.G.Evans.NE@gmail.com, on
the [MEXNC, SNCTOOLS, and the NetCDF Toolbox](http://mexcdf.sf.net) web
page.

ncvtk {#ncvtk}
-----------------------------

[Ncvtk](http://ncvtk.sourceforge.net/) is a program for exploring
planetary data stored in a NetCDF file. The NetCDF file should loosely
follow the [CF metadata
conventions](http://www.cgd.ucar.edu/cms/eaton/cf-metadata/).

Ncvtk was designed from the ground up with the aim of offering a high
degree of interactivity to scientists who have a need to explore
structured, three-dimensional, time-dependent climate data on the
sphere. A graphical user interface allows users to interact with their
data via color/transparency/contour/vector plots, apply vertical slices,
probe data, apply an external sun light, overlay hydrographic and
geopolitical data, rotate, zoom, etc. with minimal fuss.

Ncvtk is written in python and is based on the [Visualization Toolkit
(VTK)](http://public.kitware.com/VTK/). Like python and VTK, Ncvtk is
highly portable and known to run on Windows and Linux (i386, ia64,
EMT64) platforms. More information about Ncvtk is available at
<http://ncvtk.sourceforge.net>.

Ivan Shmakov's netcdf tools {#netcdf_tools}
----------------------------------------------------------

The NetCDF tools is a free software package consisting of a few tools
operating on NetCDF and, by utilizing the compatibility API, HDF4 files,
which are intended to be usable from Shell scripts.

The currently packaged tools are:

-   a couple of simple shell wrappers over the respective NetCDF
    functions (ncattget and ncattput);
-   a more sophisticated ncget tool.

The ncget tool implements functionalilty that is similar to hdp dumpsds
(for NetCDF, which lacks such a tool), or complements it in the case of
HDF4. It can be seen as a complement to the ncdump tool (included in
both the NetCDF and HDF4 distributions) as well.

This tool allows a selected part of a NetCDF variable or an HDF4
scientific data set (SDS) to be extracted in either an ASCII or binary
form, applying the transformation specified by the usual scale\_factor
and add\_offset attributes. It allows one to feed the data contained in
NetCDF variables (or HDF4 SDS) to the tools designed to operate on
either ASCII (text) or raw (binary) data.

This version of the package is the first one to be announced to the
public. It has some known bugs and limitations, but it's proved to be
quite usable. A [project
page](http://freshmeat.net/projects/netcdf-tools) on freshmeat.net. The
[source](http://waterlily.siamics.net/~ivan/src/netcdf-tools-0.1-rc1.tar.gz)
is also available.

netcdf4excel (add-in for MS Excel) {#netcdf4excel}
-----------------------------------------------------------------

Alexander Bruhns has developed [a netCDF add-in written in Visual Basic
for MS Excel](http://code.google.com/p/netcdf4excel/). This add-in
simplifies the use of NetCDF data in Excel, providing a ready to use
solution for manipulating this type of data.

For developers, the open-source (GPL V3 license) can be downloaded
directly or checked out with Mercurial.

The add-in is written in VBA 6.0 (so it won't work with Office 2010 64
bits) and is designed for Excel 2007 running with the Microsoft Windows
operating system. It supports opening netCDF classic format data with
Excel for read or write access.

More details are available on the [netcdf4excel web
site](http://code.google.com/p/netcdf4excel/).

NetCDF95 alternative Fortran API {#netcdf95}
-----------------------------------------------------------

Lionel Guez has developed and made feely available
[NetCDF95](http://web.lmd.jussieu.fr/~lglmd/NetCDF95), a new alternative
Fortran interface to the NetCDF library. Compared to the
Unidata-provided Fortran 90 netCDF interface, the NetCDF95 interface is
meant to be easier to use and more secure.

Objective-C API {#Objective-C}
---------------------------------------------

Tom Moore has an Objective-C API, available here:
[www.paleoterra.com/software](http://www.paleoterra.com/software). The
netCDF Framework is an open source (Argonne Open Source License) MacOSX
application framework that provides an Objective-C interface to the NCAR
netCDF library version 3. The framework is available both as source code
and universal compiles (works on both PPC and Intel macs). The source
code has also been compiled by users for the GNUStep environment.
Version 2 of the framework will provide classes for accessing multiple
netCDF files, working with in-memory data slabs using standard notation,
and some support for multithreading.

### Mark Tracy's Objective-C API

Mark Tracy has written [NetcdfStep](http://www.mt-se.com/nc_1.html), an
Objective-C API for netCDF that uses Objective-C Foundation Classes.

NetcdfStep is framework for using the netCDF library in object-oriented
programming with Objective-C. It now supports the full functionality of
netCDF 3.6.2.

A [complete Mac OS X
distribution](http://www.mt-se.com/pub/NetcdfStep-1.0.2.zip) including
pre-built static library and [online
documentation](http://www.mt-se.com/netcdfstep_doc/) are available.
Applications linked to this framework have no external dependencies
(other than Mac OS X itself). A [source-code only
distribution](http://www.mt-se.com/pub/NetcdfStep-GNUstep-0.6.1.tar.gz)
synced up to version 0.6.1 is available for GNUstep for use on Linux and
other Unix platforms.

Octave interface {#NCMEX}
----------------------------------------

The ARM Program has contributed NCMEX for Octave, a port of Chuck
Denham's MATLAB NCMEX to [Octave](http://www.octave.org). The calling
syntax is identical, so scripts using NCMEX in MATLAB should in theory
be portable to Octave. In order to build NCMEX, a compiled C NetCDF
library must already be installed.

In addition to the base NetCDF library interface, this package includes
a simple toolbox to automate the reading and writing of NetCDf files
within Octave using NCMEX. These tools as well as the source for NCMEX
are available from
<http://engineering.arm.gov/~sbeus/octavex/octavex.tar> (NOTE: this .tar
file contains other Octave extension functions besides NCMEX.)

Also see [Octcdf](http://ocgmod1.marine.usf.edu/octcdf/), a netCDF
toolbox for Octave.

For installation instructions, see the README file inside the .tar file.

Octave interface (Barth) {#Octave}
-------------------------------------------------

Alexander Barth has contributed the following:

Octcdf is a netCDF toolbox for [Octave](http://www.octave.org/) which
uses the same operator syntax as the [matlab netCDF
toolbox](http://mexcdf.sourceforge.net/netcdf_toolbox.html) of Charles
R. Denham. NetCDF dimensions, attributes and variables are Octave
objects and can be accessed, sliced and changed just as regular
variables. Unlike most netCDF toolboxes for matlab, it does not depend
on the NCMEX wrapper around the netCDF interface. This octave toolbox is
written in C++ calling directly the netCDF library. The octcdf toolbox
can also be used to download data from an OpenDAP server. The octcdf
source code is available at
<http://modb.oce.ulg.ac.be/mediawiki/index.php/NetCDF_toolbox_for_Octave>.
It was also included in the Octave Repository
[octave-forge](http://octave.sourceforge.net/).

OPeNDAP (formerly DODS) {#OPeNDAP}
-------------------------------------------------

The [OPeNDAP](http://opendap.org/) (formerly known as DODS) is an
Open-source Project for a Network Data Access Protocol that makes local
data and subsets of local data accessible to remote locations
independent of the local storage format. OPeNDAP also provides tools for
transforming existing applications into OPeNDAP clients, enabling them
to remotely access OPeNDAP served data. OPeNDAP is based on existing
data access tools; rather than developing a self contained system, it
makes extensive use of existing data access APIs.

OPeNDAP can be used to make netCDF data files available over the
Internet and it can also be used to adapt existing software which use
the netCDF API (by re-linking) to read data served by an OPeNDAP data
server. In principle, any program written using netCDF can be adapted to
read data from an OPeNDAP server - in other words any program which uses
netCDF can become a client in the OPeNDAP client-server system. Included
in the source and binary distributions are two freely available programs
that have already been modified (re-linked).

With a client program accessing data from a netCDF server, it is
possible to access a small subset of a large dataset over the Internet
without copying the entire dataset (as you would have to do with FTP or
AFS). The client can see changes to the netCDF dataset, e.g. when new
records are added (which would not be possible with FTP). Finally, the
client can also access cross-sections of variable data without paging
large amounts of data across the network (as you would have to do with
NFS, for example).

OPeNDAP software is freely available in both source form or binary form
for selected platforms.

OpenDX {#OpenDX}
-------------------------------

[OpenDX](http://www.opendx.org/about.html) (formerly IBM Data Explorer,
also known as simply DX) is a general-purpose software package for data
visualization and analysis. It employs a data-flow driven client-server
execution model and provides a graphical program editor that allows the
user to create a visualization using a point and click interface.
DX runs on 7 major UNIX platforms as well as Windows 95/NT and is
designed to take full advantage of multi-processor systems from IBM, SGI
and Sun.

DX is built upon an internal data model, which describes and provides
uniform access services for any data brought into, generated by, or
exported from the software. This data model supports a number of
different classes of scientific data, which can be described by their
shape (size and number of dimensions), rank (e.g., scalar, vector,
tensor), type (float, integer, byte, etc. or real, complex, quaternion),
where the data are located in space (positions), how the locations are
related to each other (connections), aggregates or groups (e.g.,
hierarchies, series, composites, multizone grids, etc.). It also
supports those entities required for graphics and imaging operations
within the context of Data Explorer. Regular and irregular, deformed or
curvilinear, structured and unstructured data as well as "missing" or
invalid data are supported.

The details of the data model are hidden at the user level. As a result
DX operations or modules are polymorphic and appear typeless. The DX
Import module, which reads data for use within Data Explorer directly
utilizes data in netCDF as well as other formats (e.g., HDF, CDF). One
or more variables may be selected as well as step(s) of a time series.
Data in conventional netCDFs are directly imported. Since the DX data
model is more comprehensive than the netCDF data model, a methodology to
extend netCDF via attribute conventions (e.g., for unstructured meshes,
non-scalar data and hierarchies) for use with Data Explorer is
available.

DX supports a number of realization techniques for generating renderable
geometry from data. These include color and opacity mapping (e.g., for
surface and volume rendering), contours and isosurfaces, histograms,
two-dimensional and three-dimensional plotting, surface deformation,
etc. for scalar data. For vector data, arrow plots, streamlines,
streaklines, etc. are provided. Realizations may be annotated with
ribbons, tubes, axes, glyphs, text and display of data locations, meshes
and boundaries. Data probing, picking, arbitrary surface and volume
sampling, and arbitrary cutting/mapping planes are supported.

DX supports a number of non-graphical functions such as point-wise
mathematical expressions (e.g., arithmetic, transcendental, boolean,
type conversion, etc.), univariate statistics and image processing
(e.g., transformation, filter, warp, edge detection, convolution,
equalization, blending, morphological operations, etc.). Field/vector
operations such as divergence, gradient and curl, dot and cross
products, etc. are provided. Non-gridded or scattered data may be
interpolated to an arbitrary grid or triangulated, depending on the
analysis requirements. The length, area or volume of various geometries
may also be computed. Tools for data manipulation such as removal of
data points, subsetting by position, sub/supersampling, grid
construction, mapping, interpolation, regridding, transposition, etc.
are available.

Tools for doing cartographic projections and registration as well as
earth, space and environmental sciences examples are available at
Cornell University via info.tc.cornell.edu. Also see the [ncdx](#ncdx)
tool for making netCDF files OpenDX compliant.

Panoply {#Panoply}
---------------------------------

[Panoply](http://www.giss.nasa.gov/tools/panoply/) is an application
that plots geo-gridded and other arrays from netCDF, HDF, GRIB, and
other datasets. Features include:

-   Slice and plot geo-gridded latitude-longitude, latitude-vertical,
    longitude-vertical, or time-latitude arrays from larger
    multidimensional variables.
-   Two arrays may be combined in one plot by differencing, summing, or
    averaging.
-   Lon-lat data may be plotted as global maps (using any of over 75 map
    projections) or as zonal average plots.
-   Overlay continent outlines or masks on lon-lat plots.
-   Use your favorite CPT, GGR, PAL, or ACT color table for scale
    colorbar.
-   Save plots to disk in GIF, JPEG, PNG or TIFF bitmap images or as PDF
    or PostScript graphics files.
-   Export lon-lat map plots in KMZ format.
-   Export animations as AVI or MOV video or as a collection of
    invididual frame images.
-   Explore remote THREDDS and OpenDAP catalogs and open datasets served
    from them.

Panoply requires that your computer have a Java SE 6 runtime
environment, or better, installed.

Panoply is developed at the NASA Goddard Institute for Space Studies.
Questions and suggestions should be directed to [Dr. Robert B.
Schmunk](http://www.giss.nasa.gov/staff/rschmunk.html).

Parallel-NetCDF {#Parallel-NetCDF}
-------------------------------------------------

A group of researchers at Northwestern University and Argonne National
Laboratory (Jianwei Li, Wei-keng Liao, Alok Choudhary, Robert Ross,
Rajeev Thakur, William Gropp, and Rob Latham) have designed and
implemented a new [parallel interface for writing and reading netCDF
data](http://www.mcs.anl.gov/parallel-netcdf/), tailored for use on high
performance platforms with parallel I/O. The implementation builds on
the MPI-IO interface, providing portability to most platforms in use and
allowing users to leverage the many optimizations built into MPI-IO
implementations. Testing so far has been on Linux platforms with ROMIO
and IBM SP machines using IBM's MPI.

Documentation and code for Parallel-NetCDF is now available for testing.
Although a few interfaces are not implemented yet, the current
implementation is complete enough to provide significant I/O performance
improvements on parallel platforms, as described in a [technical
report](ftp://info.mcs.anl.gov/pub/tech_reports/reports/P1048.pdf).
Users are invited to test Parallel-NetCDF in their applications.

Paraview and vtkCSCSNetCDF {#Paraview}
-----------------------------------------------------

<http://www.paraview.org/>

ParaView is an application designed with the need to visualize large
data sets in mind. The goals of the ParaView project include the
following:

-   Develop an open-source, multi-platform visualization application.
-   Support distributed computation models to process large data sets.
-   Create an open, flexible, and intuitive user interface.
-   Develop an extensible architecture based on open standards.

ParaView runs on distributed and shared memory parallel as well as
single processor systems and has been successfully tested on Windows,
Linux and various Unix workstations and clusters. Under the hood,
ParaView uses the Visualization Toolkit as the data processing and
rendering engine and has a user interface written using a unique blend
of Tcl/Tk and C++.

A vtk/ParaView reader for netCDF files can be found here.

Perl interfaces {#Perl}
--------------------------------------

There are two netCDF interfaces for Perl:
-   [PDL::NetCDF](http://search.cpan.org/~dhunt/PDL-NetCDF-4.05/netcdf.pd),
    Doug Hunt's perl interface which uses the PDL (perl data language)
    extension.
-   [NetCDFPerl](/software/netcdf-perl/), Steve Emmerson's extension
    module, based on version 2 of the netCDF package. Uses perl lists
    for representing netCDF variables.

PolyPaint+ {#PolyPaint}
---------------------------------------

[PolyPaint+](http://lasp.colorado.edu/polypaint/home.html) is an
interactive scientific visualization tool that displays complex
structures within three-dimensional data fields. It provides both color
shaded-surface display and simple volumetric rendering in either index
or true color. For shaded surface rendering, the PolyPaint+ routines
first compute the polygon set that describes a desired surface within
the 3D data volume. These polygons are then rendered as continuously
shaded surfaces. PolyPaint+ contains a wide variety of options that
control lighting, viewing, and shading. Objects rendered volumetrically
may be viewed along with shaded surfaces. Additional data sets can be
overlaid on shaded surfaces by color coding the data according to a
specified color ramp. 3D visualizations can be viewed in stereo for
added depth perspective.
Currently supported 3D visualizations are the following:

-   Shaded isosurface
-   Transparent contour shells or isosurfaces at varying levels
-   Volumetric or density plot
-   Planes
-   Contour ribbons
-   Topographic surface from 2D geographic data sets

3D data volumes may be sliced in the X, Y, or Z plane using an
interactive cutting plane. A cross section of the data volume can be
viewed in a 2D window as a 2D contour plot, a vector plot, a raster
image or a combination of these options superimposed. Map outlines can
be used as a background for 2D cross section plots of geographic data.
All data is projected according to the coordinates specified by the user
for the cross section window.

The user interface provides direct manipulation tools for specifying the
eye position, center of view, light sources, and color ramps. Subsetting
of data can be done easily by selecting the data by index or geographic
coordinate. On-line contextual help provides easy access to more detail
about the software. Tutorials which range from very simple
visualizations to complex combinations of data sets provide the user
with a quick learning tool.

Currently PolyPaint+ accepts only data which is in the NetCDF file
format. A file conversion utility which converts from raw binary data to
netCDf is a part of the application.

PolyPaint+ is a joint effort of the University of Colorado and NCAR
(National Center for Atmospheric Research) funded by the NASA AISRP
program. A beta version of PolyPaint+ is currently available free of
charge using FTP or for a nominal fee which would cover tape
distribution. A license agreement must be signed in order to use it.

You may order by...

-   TELEPHONE : 303-492-7289 (Margi Klemp) : 303-497-8159 (Bill Boyd)
-   U.S. MAIL :

                Margi Klemp
                University of Colorado / LASP
                1234 Innovation Dr.
                Boulder, CO 80303
                USA

-   E-MAIL : margi@aries.colorado.edu

Pomegranate {#Pomegranate}
-----------------------------------------

The P9E Team at NASA JPL has developed
[Pomegranate](http://pomegranate.jpl.nasa.gov/), a python application
that "webifies" science data files. Supported formats include netCDF,
HDF4, HDF5, GRIB and FITS.

Pomegranate can be installed on web servers as either a WSGI or CGI
application to provide webification (w10n) services. To learn more about
w10n of science data files, please visit <http://webification.org/>. A
brief [help](http://pomegranate.jpl.nasa.gov/test/help.txt) document
describes how to use the [demo
directory](http://pomegranate.jpl.nasa.gov/test) to browse or download
metadata or data in netCDF, JSON, or other formats by clicking on data
folder and document icons.

Pomegranate can also be used as a standalone library or command line
application. This greatly simplifies the retrieval of metadata and data
from files in supported formats.

Pomegranate is open source software and can be downloaded from
<http://www.openchannelsoftware.com/projects/Pomegranate/>.

PyNGL and PyNIO {#PyNGL}
---------------------------------------

NCAR's Computational and Information Systems Laboratory has developed
[PyNGL](http://www.pyngl.ucar.edu/), a python package for scientific
visualization and data analysis and
[PyNIO](http://www.pyngl.ucar.edu/Nio.shtml), a Python package
supporting access to a variety of data formats using an interface
modelled on netCDF.

Python interfaces {#Python}
------------------------------------------

Python is an interpreted, object-oriented language that is supported on
a wide range of hardware and operating systems. Python information and
sources can be obtained from <http://www.python.org/>. There are now
several netCDF interfaces for Python.

Jeff Whitaker of the NOAA Earth System Research Lab has developed a
netCDF-4 module for python: <http://code.google.com/p/netcdf4-python/>.
Most new features of netCDF-4 are implemented, such as multiple
unlimited dimensions, groups and zlib data compression. All the new
numeric data types (such as 64-bit and unsigned integer types) are
implemented. Compound and variable length (vlen) data types are
supported, but the enum and opaque data types are not. Mixtures of
compound and vlen data types (compound types containing vlens, and vlens
containing compound types) are not supported.

[xray](#xray) is a higher-level interface that uses netcdf4-python
internally to implement a pandas-like package for N-D labelled arrays
for scientific data.

André Gosselin of the Institut Maurice-Lamontagne, Péches & Océans
Canada, has implemented pycdf, a new Python interface to the netCDF
library. It is available from <http://pysclint.sourceforge.net/pycdf/>,
where you will find the install files, installation instructions,
extensive documentation in text and html format, and examples. pycdf
requires the Numeric python package, and installs through the simple
"python setyp.py install" command.

Bill Noon (noon@snow.cit.cornell.edu) has implemented another netCDF
Python module that allows easy creation, access, and browsing of netCDF
data. The bindings also use the [udunits library](/software/udunits/) to
do unit conversions. More information and source for Noon's Python
netCDF module are available from
<http://snow.cit.cornell.edu/noon/ncmodule.html>.

The package from Konrad Hinsen has been integrated into his
[ScientificPython](https://sourcesup.cru.fr/projects/scientific-py/)
package.

Dave Brown of NCAR's Computational and Information Systems Laboratory
has developed [PyNIO](http://www.pyngl.ucar.edu/Nio.shtml), a Python
package that allows read and/or write access to a variety of data
formats using an interface modelled on netCDF. Currently supported
formats include netCDF, HDF4, GRIB1 and GRIB2 (read only), and HDF-EOS 2
Grid and Swath data (read only).

Vicente Galiano of Miguel Hernandez University has developed a Python
interface to PnetCDF. This Python's package called "PyPnetCDF" allows
access to NetCDF files using MPI and the library pnetCDF developed by
http://www.mcs.anl.gov/parallel-netcdf/. The tools are very similar to
Konrad Hinsen's NetCDF package to Python but can read and write in a
parallel way. For more information, see:
<http://www.pyacts.org/pypnetcdf>.

Pupynere (PUre PYthon NEtcdf REader) Roberto
De Almeida has developed
[pupynere](http://pypi.python.org/pypi/pupynere/), a PUre PYthon NEtcdf
REader that allows read-access to netCDF files using the same syntax as
the Scientific.IO.NetCDF Python module. Even though it's written in
Python, the module is up to 40% faster than Scientific.IO.NetCDF and
pynetcdf.

R interface {#R}
-------------------------------

The R Project for Statistical Computing has developed
[R](http://www.R-project.org/), a language and environment for
statistical computing and graphics. It provides a wide variety of
statistical and graphical techniques, including linear and nonlinear
modelling, statistical tests, time series analysis, classification, and
clustering.

David Pierce has contributed the [ncdf4
package](http://cran.r-project.org/web/packages/ncdf4/index.html) for
reading netCDF data into R and for creating new netCDF dimensions,
variables, and files, or manipulating existing netCDF files from R.

Pavel Michna has contributed another package,
[RNetCDF](http://cran.r-project.org/web/packages/RNetCDF/index.html),
that also provides access to netCDF data and to udunits calendar
functions from R.

Robert Hijmans (with additional contributors) has created the [R raster
package](http://cran.r-project.org/web/packages/raster/index.html) for
geographic data analysis and modeling. The raster package can be used
for reading, writing, manipulating, analyzing and modeling gridded
spatial data. The package is especially useful for large datasets that
don't fit into memory, because data is processed in chunks. See
[Introduction to the 'raster'
package](http://cran.r-project.org/web/packages/raster/vignettes/Raster.pdf)
for more information.

Quantum GIS (QGIS) {#QGIS}
-----------------------------------------

[Quantum GIS](http://www.qgis.org/) (QGIS) is an Open Source Geographic Information System (GIS) licensed under the GNU General Public License. QGIS is an official project of the Open Source Geospatial Foundation (OSGeo). It runs on Linux, Unix, Mac OSX, and Windows and supports numerous vector, raster, and database formats and functionalities. QGIS supports a desktop, browser, server, and client for viewing, editing, analysis, serving, and accessing data. Its server complies with the OGC WMS 1.3 standard. In addition to PostGIS and SpatiaLite formats, it can access data in vector formats supported by the OGR library as well as most raster formats supported by the GDAL library, including netCDF. For a more detailed list of features of the QGIS desktop, browser, server, and client, see the [QGIS features page](http://www.qgis.org/en/about-qgis/features.html).

QGIS also supports displaying netCDF using the [Crayfish plugin in QGIS](http://www.lutraconsulting.co.uk/blog/2016/04/14/crayfish-2-2/).  The project repository may be found [here](https://github.com/lutraconsulting/qgis-crayfish-plugin).

Ruby interface {#Ruby}
-------------------------------------

A group at the Research Institute for Sustainable Humanosphere (RISH) of
Kyoto University has developed a [netCDF interface for
Ruby](http://www.gfd-dennou.org/arch/ruby/products/ruby-netcdf/), an
interpreted, object-oriented scripting language. This interface is
intended to cover all the functionality of the C library for netCDF.
Also available are combination functions such as iterators (which offer
abstract ways to scan files and variables). Numeric arrays are handled
by the "NArray" multi-dimensional array class, which is becoming the de
facto standard multi-dimensional array for Ruby. See also the Ruby-based
[GPhys software and Gfdnavi tool](#Gfdnavi) for accessing GRIB, GrADS,
and netCDF data uniformly.

More information about Ruby is available from the [Ruby web
site](http://www.ruby-lang.org/).

Scientific DataSet (SDS) Library {#SDS}
------------------------------------------------------

The [Scientific DataSet Library and Tools
project](http://sds.codeplex.com), developed jointly by Microsoft
Research Cambridge and Moscow State University, is aimed at manipulation
and visualization of multidimensional data sets.

Scientific DataSet (or SDS in short) is a .NET class library for
manipulating scientific data and their metadata. SDS provides a unified
API for convenient access to various data storages. Three types of
storages are supported by the first release: NetCDF files, CSV text
files and volatile in-memory datasets. SDS uses native NetCDF library
built from version 4.0.1 both for 32 and 64-bit Windows platforms. New
storage types can be added to SDS infractructure as plugins. Support for
accessing TIFF image files from SDS as 2D arrays will be available soon
as a separate CodePlex project.

Three applications are built on top of SDS:

-   sds command line utility. It allows users to examine data set
    schema, copy data sets, modify their metadata.
-   DataSetViewer application for visualization of data sets.
    DataSetViewer is both a standalone application and Windows
    Presentation Foundation Control that can be built into your
    applications. DataSetViewer has support for interactive slicing of
    multidimensional data along any dimension.
-   DataSetEditor add-in for Microsoft Office Excel. DataSetEditor
    provides ability to view and modify the contents of any data set as
    Excel worksheets.

You can read the Getting Started document at
<http://sds.codeplex.com/Project/Download/FileDownload.aspx?DownloadId=127282>
for a more detailed introduction to the Scientific DataSet software. A
Windows Installation package for SDS binaries along with DataSet Viewer
and DataSet Editor are available also. You can also build core class
libraries and the sds utility under Mono. You may use, copy, and
reproduce this software for any non-commercial purpose. For further
details see license at <http://sds.codeplex.com/license>.

The SDS project is in beta phase and keeps evolving. You are welcome to
join discussions or report issues at the CodePlex site:
<http://sds.codeplex.com>.

sciNetCDF {#scinetcdf}
-------------------------------------------------------------
[sciNetCDF](https://atoms.scilab.org/toolboxes/scinetcdf)

In the context of the IASI-NG project, CNES is responsible for the development
of a Scilab/NetCDF4 interface, which CNES wanted to make available to the entire
scientific community.

The toolbox sciNetCDF is the result of this collaboration. It can read and write
NetCDF files of any version (version 4 of the format is used by default for
writing).

The toolbox provides high level functions to read/write NetCDF files natively in
Scilab in a friendly manner (data is converted automatically from Scilab to
NetCDF and inversely).
These functions are:
- nccreate
- ncwrite
- ncread
- ncwriteatt
- ncreadatt
- ncdisp

It provides also a low level interface to all the NetCDF C library functions


Apache Spatial Information System (SIS) {#SIS}
-------------------------------------------------------------

[Apache Spatial Information System
(SIS)](https://builds.apache.org/job/sis-trunk/site/index.html) is a
Java library for developing geospatial applications. SIS enables
representation of coordinates for searching, data clustering, archiving,
or any other relevant spatial needs. The library is an implementation of
GeoAPI 3.0 interfaces and can be used for desktop or server
applications.

SIS provides data structures for geographic data and associated metadata
along with methods to manipulate those data structures. The SIS metadata
module forms the base of the library and enables the creation of
metadata objects which comply with the ISO 19115 metadata model and
which can be read from or written to ISO 19139 compliant XML documents.
The SIS referencing module will enable the construction of geodetic data
structures for geospatial referencing based on the ISO 19111 model such
as axis, projection and coordinate reference system definitions, along
with the associated operations which enable the mathematical conversion
of coordinates between different systems of reference. The SIS storage
modules will provide a common approach to the reading and writing of
grid coverages applicable to simple imagery and multidimensional data
structures.

SIS supports creating ISO 19115 metadata from metadata in a netCDF store
from a given file, URL, stream, or NetcdfFile object. SIS netCDF storage
is intended to be a bridge between NetCDF Climate and Forecast (CF)
conventions and ISO 19115 metadata.

SIS is under development as an Apache project. Release 0.3 is currently
available for download.

Tcl/Tk interfaces {#TclTk}
------------------------------------------

Dan Schmitt has developed [cdftcl](http://cnrit.tamu.edu/rsg/cdftcl/), a
[Tcl/Tk](http://www.scriptics.com/) interface for netCDF. It allows the
use of "wildcards" (\*) or ranges (1-4) in the subscript notation, and
use of name references instead of variable IDs. Contact dan@computer.org
for more information.

Tcl-nap {#Tcl-nap}
---------------------------------

[Tcl-nap](http://tcl-nap.sourceforge.net) (n-dimensional array
processor) is a loadable extension of Tcl which provides a powerful and
efficient facility for processing data in the form of n-dimensional
arrays. It has been designed to provide an array-processing facility
with much of the functionality of languages such as
[APL](http://www.acm.org/sigapl/), Fortran-90, [IDL](#IDL),
[J](http://www.jsoftware.com/), [matlab](http://www.mathworks.com), and
[octave](http://www.octave.org/).

Support is provided for data based on n-dimensional grids, where the
dimensions correspond to continuous spatial coordinates. There are
interfaces to the HDF and netCDF file formats commonly used for such
data, especially in Earth sciences such as Oceanography and Meteorology.

The internal data structure is called a NAO (n-dimensional array object)
and contains similar information to that of HDF SDSs and netCDF
variables.

Tcl-nap was developed as part of the [CSIRO CAPS
project](http://www.dar.csiro.au/rs/avhrr_processing_software.htm), but
can be loaded and used without the (satellite oriented) CAPS extension.

Visual Basic and VB.net interfaces {#VB}
-------------------------------------------------------

Carsten Wieczorrek has developed code in VB 6 to export chromatographic
data into the netcdf/ANDI format. The application writes netCDF files
that can be read by CHROMELEON, for example. For others interested in
programming with netcdf.dll from VB 6, see Wieczorrek's web page on
[netCDF and VB 6.0](http://www.mn-net.com/netcdf_vb6) and for VB.net,
see [netCDF and VB.net](http://www.mn-net.com/netcdf_vbnet).

VisAD {#VisAD}
-----------------------------

[VisAD](http://www.ssec.wisc.edu/~billh/visad.html) is a Java class
library for interactive and collaborative visualization and analysis of
numerical data. It combines:
-   The use of pure Java for platform independence and to support data
    sharing and real-time collaboration among geographically distributed
    users. Support for distributed computing is integrated at the lowest
    levels of the system using Java RMI distributed objects.
-   A general mathematical data model that can be adapted to virtually
    any numerical data, that supports data sharing among different
    users, different data sources and different scientific disciplines,
    and that provides transparent access to data independent of storage
    format and location (i.e., memory, disk or remote). The data model
    has been adapted to netCDF, FITS, HDF-EOS, McIDAS, Vis5D, GIF and
    JPEG file formats.
-   A general display model that supports interactive 3-D, data fusion,
    multiple data views, direct manipulation, collaboration, and virtual
    reality. The display model has been adapted to Java3D and Java2D and
    used in an ImmersaDesk virtual reality display.
-   Data analysis and computation integrated with visualization to
    support computational steering and other complex interaction modes.
-   Support for two distinct communities: developers who create domain-
    specific systems based on VisAD, and users of those domain-specific
    systems. VisAD is designed to support a wide variety of user
    interfaces, ranging from simple data browser applets to complex
    applications that allow groups of scientists to collaboratively
    develop data analysis algorithms.
-   Developer extensibility in as many ways as possible.

VisAD was written by programmers at the [SSEC Visualization
Project](http://www.ssec.wisc.edu/~billh/vis.html) at the University of
Wisconsin-Madison [Space Science and Engineering
Center](http://www.ssec.wisc.edu/), and the [Unidata Program
Center](/index.html).

WebWinds {#WebWinds}
-----------------------------------

[WebWinds](http://www.openchannelsoftware.com/projects/WebWinds/) is a
free Java-based science visualization and analysis package. In addition
to several new analysis tools, the current fourth version does automatic
scripting. This allows

1.  a user to rapidly and automatically create and store a session,
    either for his own use, or for use by a collaborator on another
    machine;
2.  a data provider to automatically create a specialized analysis
    environment which can be downloaded (as a small script file) along
    with a dataset from a Website; and
3.  realtime collaboration or sharing of sessions over (even
    low-bandwidth) networks, including the Internet.

This scripting requires no knowledge of the scripting language syntax.
Several sample script files are included with the distribution.

In addition, this version contains a capability to geo-reference some
data and to read ASCII data in tabular format. Also new is the ability
to output data in numerical form (e.g. NetCDF) and a context sensitive,
integrated help system.

As with earlier versions, data in several different formats, including
NetCDF, can be read in easily from your local machine or from the Web.
In addition, most data can be subset or subsampled on load, making it
possible to visualize very large multidimensional and/or multispectral
datasets. The package includes several step-by-step examples.
Installation of the software (including Java) on the PC or Mac is a
process requiring one file to be downloaded and opened. If you need help
getting started, a remote tutorial is available once you've downloaded
the package.

WebWinds is \`point and click' rather than language driven and it runs
well on Unix, Windows (95/98/NT) and Mac platforms. It currently
requires JDK 1.1. To download a copy of this release, go to
<http://www.sci-conservices.com/rel4/webpage/wwhome.html>

xray (Python N-D labelled arrays) {#xray}
--------------------------------------------------------

[xray](http://xray.readthedocs.org/en/stable/index.html) is an open
source project and Python package that aims to bring the labeled data
power of [pandas](http://pandas.pydata.org/) to the physical sciences,
by providing N-dimensional variants of the core pandas data structures,
Series and DataFrame: the xray DataArray and Dataset.

xray adopts the [Common Data
Model](http://www.unidata.ucar.edu/software/thredds/current/netcdf-java/CDM)
for self-describing scientific data in widespread use in the Earth
sciences (e.g., netCDF and OPeNDAP): xray.Dataset is an in-memory
representation of a netCDF file.

xray is being developed by Stephan Hoyer, Alex Kleeman, and [other
contributors](https://github.com/xray/xray/graphs/contributors).

Zebra {#Zebra}
-----------------------------

[Zebra](http://www.atd.ucar.edu/rdp/zebra.html) (formerly named Zeb) is
a system for data ingest, storage, integration and display, designed to
operate in both real time and postprocessing modes. Zebra was developed
by Jonathan Corbet and others in NCAR's [Research Data
Program](http://www.atd.ucar.edu/rdp/rdp_home.html).
Zebra's primary use is for the superpositioning of observational data
sets (such as those collected by satellite, radar, mesonet and aircraft)
and analysis products (such as model results, dual-Doppler synthesis or
algorithm output). Data may be overlaid on a variety of display types,
including constant altitude planes, vertical cross-sections, X-Y graphs,
Skew-T plots and time-height profiles. The fields for display, color
tables, contour intervals and various other display options are defined
using an icon based user-interface. This highly flexible system allows
scientific investigators to interactively superimpose and highlight
diverse data sets; thus aiding data interpretation.

Data handling capabilities permit external analysis programs to be
easily linked with display and data storage processes. The data store
accepts incoming data, stores it on disk, and makes it available to
processes which need it. An application library is available for data
handling. The library functions allow data storage, retrieval and
queries using a single applications interface, regardless of the data's
source and organization. NetCDF data that conforms to Zebra conventions
is supported by this interface.

Zebra is currently available to the university research community
through the NCAR/ATD Research Data Program. Email requests to
rdp-support@atd.ucar.edu. More information is on the web page
http://www.atd.ucar.edu/rdp/zebra.html.

------------------------------------------------------------------------

User-Contributed Software {#user}
================================================

Unidata makes available a separate
[catalog](/software/netcdf/Contrib.html) to a
[directory](ftp://ftp.unidata.ucar.edu/pub/netcdf/contrib/) of freely
available, user-contributed software and documentation related to the
netCDF library. This software may be retrieved by anonymous FTP. We
haven't necessarily used or tested this software; we make it available
"as is".

The criteria for inclusion in the netcdf/contrib/ directory of
user-contributed software are:

-   General usefulness to a significant part of the netCDF community
-   Small size
-   Infrequent need for updates
-   Free availability

------------------------------------------------------------------------

Commercial or Licensed Packages {#commercial}
===============================

ASA ViewNcDap {#ViewNcDap}
-----------------------------------------

Applied Science Associates, Inc. has made the ASA View NC/Dap
application freely available for
[download](http://www.asascience.com/downloads). ViewNcDap is a
stand-alone research-based tool (with included demonstration data) that
allows a user to visualize four dimensional NetCDF and OPeNDAP data.
ViewNcDap is a Windows application that includes temporal/time step
functionality for viewing animations of data that include temporal
information. The application may be used to visualize a variety of
time-varying geospatial scientific data in a simple map framework. It
handles CF conventions and includes some aliasing features that could
permit additional formats to be read. It should not be considered a GIS
system, but is used to quickly preview a variety of data on a simple
map. Data may also be filtered and saved to a local netCDF file.

Avizo {#Avizo}
-----------------------------

[Avizo](http://www.avizo3d.com/) software is a powerful tool for 3D data
visualization and analysis. It offers a comprehensive feature set that
addresses visualization, processing, analysis, communication and
presentation. [Avizo Green
Edition](http://www.vsg3d.com/vsg_prod_avizo_green.php) includes an
advanced set of features dedicated to climate, oceanography,
environmental or earth-mapped data. It provides high-level support for
the netCDF format, a dedicated Earth visualization module, and a set of
advanced geographical projections applicable to a wide range of fast 2D
and 3D data representations.

For more information, see [www.avizo3d.com](http://www.avizo3d.com/).

AVS {#AVS}
-------------------------

[AVS](ftp://testavs.ncsc.org/avs/Info/WHAT_IS_AVS) (Application
Visualization System) is a visualization application software and
development environment. An AVS module has been written that allows
multi-dimensional netCDF data sets to read into AVS as uniform or
rectilinear field files. The AVS user can point and click to specify the
name of the variable in the selected netCDF file, as well as selecting
the hyperslab. If 1D coordinate variables exist (a variable that has the
same name as a dimension) then the coordinate variable will be used to
specify the coordinates of resulting rectilinear field file. If no
coordinate variable exists, then the resulting field file will be
uniform. Once in AVS, there are hundreds of analysis and display modules
available for image processing, isosurface rendering, arbitrary slicing,
alpha blending, streamline and vorticity calculation, particle
advection, etc. AVS runs on many different platforms (Stardent, DEC,
Cray, Convex, E and S, SET, Sun, IBM, SGI, HP, FPS and WaveTracer), and
it has a flexible data model capable of handling multidimensional data
on non-Cartesian grids.
The module source code and documentation is available from the
[International AVS Center](http://iac.ncsc.org/), in the
<ftp://testavs.ncsc.org/avs/AVS5/Module_Src/data_input/read_netcdf/>
directory.

See also the information on [DDI](#DDI) for another way to use netCDF
data with AVS.

Barrodale UFI {#BCS-UFI}
---------------------------------------

[Barrodale Computing Services Ltd.](http://www.barrodale.com) (BCS) has
developed a product that addresses one of the main objections heard from
"technologists" (e.g., scientists, engineers, and other researchers) who
avoid using databases to manage their data: "my very large data files
are too cumbersome/difficult/slow/costly to load into a database". In
addition to netCDF, these files come in a variety of formats (HDF5,
GRIB, NITFS, FITS, etc.).

This BCS product is called the [Universal File Interface
(UFI)](http://www.barrodale.com/bcs/universal-file-interface-ufi); it's
a database extension based on the IBM Informix Virtual Table Interface
(VTI). *(Please continue reading even if you don't have Informix running
on your system, because IBM has just made available, at no charge, the
[Innovator-C
Edition](http://www-01.ibm.com/software/data/informix/innovator-c-edition/)
of Informix.)* A demo that uses UFI to access wind speeds can be seen
[here](http://www.barrodale.com/bcs/universal-file-interface-animation).

VTI is a technology that supports making external datasets appear as
tables to SQL queries and statements. UFI is a BCS database extension
for delivering the contents of external data files as though they were
rows in a database table. UFI makes a file look like a set of database
tables, so "UFI managed tables" are actually virtual database tables.
Consequently, users of UFI can perform SQL queries on their files
without having to first load them into a database.

DioVISTA/Storm {#DioVISTAStorm}
-----------------------------------------------

[DioVISTA/Storm](http://www.hitachi-power-solutions.com/products/product03/p03_61.html)
is a commercial software package that visualizes content of netCDF files
as a time series of grids, isosurfaces, and arrows on a 3D virtual
earth. Its user interface is similar to standard 3D earth visualizing
software. It displays OGC KML files, Shapefiles, and online map
resources through OGC Web Tile Map Services (WTMS). It supports CF
Conventions version 1.6 (lon-lat-alt-time axis and trajectory). Its
first version was released on Aug 5 2014.

Environmental WorkBench {#Environmental_WorkBench}
-----------------------------------------------------------------

[SuperComputer Systems Engineering and Services
Company](http://www.ssesco.com/) (SSESCO) has developed the
[Environmental WorkBench](http://www.ssesco.com/files/ewb.html) (EWB),
an easy to use visualization and analysis application targeted at
environmental data. The EWB currently has numerous users in the fields
of meteorological research, air quality work, and groundwater
remediation.
EWB system features include:

-   Random access file structure using the netCDF-based public domain
    MeRAF file system with support for gridded, discrete (non-grid-based
    observation), and particle types
-   Support for geo-referenced or Cartesian coordinate systems
-   Object oriented Graphical User Interface (GUI) that is very easy to
    use
-   Tools for converting model and observational data sets and data
    writers to netCDF
-   Interactive rotation/translation of scenes in 3D space
-   Time sequencing controls to step forward/backward, animate
    sequentially, or go to a chosen time step; including multiple
    asynchronous or non-uniform time steps
-   Interactive slicers to select cross sections through 3D data sets
-   Display operators available on the slices, including
    -   Contour lines with selectable contour levels
    -   Color shading by data value with variable transparency level
    -   Arrow and streamline representation for vector quantities
    -   Positional reference lines at user selected intervals
    -   Color coded shapes at each grid node
-   Multiple 3D isosurfaces at selected parameters and values with
    variable transparency
-   Display of particle positions with coloring by type, height, and
    source
-   Display of discrete data using colored spheres and labels for scalar
    data and arrows for vectors (with arrowheads or meteorological
    style)
-   Multiple user definable color maps to which isosurface and colored
    field shading may be separately assigned
-   On screen annotation for generation of report ready figures
-   Image export in any of the common image formats (gif, tiff,
    encapsulated postscript, etc.)
-   Graceful handling of missing or bad data values by all the graphics
    rendering routines
-   Automatic data synchronization to allow automatic screen updating as
    new data arrives in real-time from a model or set of sensors
-   Two and three dimensional interpolation from scattered observations
    to a grid, using the Natural Neighbor Method. This robust volume
    based method yields results far superior to distance weighting
    schemes.

Systems currently supported include Win95, WinNT, OS/2, IBM RS/6000,
Silicon Graphics, HP and SUN workstations.

SSESCO has implemented a meta-file layer on top of the netCDF library,
called MeRAF. It handles multiple netCDF files as well as automatic
max-min calculations, time-varying gridded, particle, and discrete data,
logical groupings for discrete data, and an overall simplified and
flexible interface for storing scientific data. MeRAF is being used by
the DOE at the Hanford-Meteorological Site for observational data and
will be used for their weather-modeling.

ESRI {#ESRI}
---------------------------

[ESRI ArcGIS](http://www.esri.com/software/arcgis/index.html) version
9.2 and later support [accessing netCDF time-based and multidimensional
data](http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?TopicName=An_overview_of_data_support_in_ArcGIS)
that follows CF or COARDS conventions for associating spatial locations
with data. A selected slice of netCDF data may be displayed in ArcGIS as
a raster layer, feature layer, or table. You can also drag a netCDF file
from Windows Explorer and drop it in an ESRI application such as ArcMap.

FME {#FME}
-------------------------

[FME](http://www.safe.com/fme), developed by [Safe Software
Inc.](http://www.safe.com), is a tool for transforming data for exchange
between over [300 different formats and
models](http://www.safe.com/fme/format-search/), including netCDF. FME's
read and write support for netCDF allows users to move data into the
netCDF common standard, regardless of its source, and conversely enables
end-users to consume netCDF data for use in their preferred systems. For
more information visit <http://www.safe.com/fme>.

HDF Explorer {#HDF-Explorer}
-------------------------------------------

[HDF Explorer](http://www.space-research.org/) is a data visualization
program that reads the HDF, HDF5 and netCDF data file formats (including
netCDF classic format data). HDF Explorer runs in the Microsoft Windows
operating systems.

HDF Explorer offers a simple yet powerful interface for the
visualization of HDF and netCDF data. The data is just a click of the
mouse away. Data is first viewed in a tree-like interface, and then
optionally loaded and visualized in a variety of ways. HDF Explorer
features include fast access to data, grid, scalar and vector views. It
also allows exporting your data either as an ASCII text file or a bitmap
image.

IDL Interface {#IDL}
-----------------------------------

[IDL](http://www.exelisvis.com/ProductsServices/IDL.aspx) (Interactive
Data Language) is a scientific computing environment, developed and
supported by [Excelis Visual Information
Solutions](http://www.exelisvis.com/), that combines mathematics,
advanced data visualization, scientific graphics, and a graphical user
interface toolkit to analyze and visualize scientific data. Designed for
use by scientists and scientific application developers, IDL's
array-oriented, fourth-generation programming language allows you to
prototype and develop complete applications. IDL now supports data in
netCDF format.
As an example, here is how to read data from a netCDF variable named GP
in a file named "data/aprin.nc" into an IDL variable named gp using the
IDL language:

       id = ncdf_open('data/april.nc')
        ncdf_varget,id, ncdf_varid( id, 'GP'), gp

Now you can visualize the data in the gp variable in a large variety of
ways and use it in other computations in IDL. You can FTP a demo version
of IDL, including the netCDF interface, by following the instructions in
pub/idl/README available via anonymous FTP from gateway.rsinc.com or
boulder.colorado.edu.
Other software packages that use or interoperate with IDL to access
netCDF data includes [ARGOS](#ARGOS), [CIDS Tools](#CIDS_Tools),
[DDI](#DDI), [HIPHOP](#HIPHOP), [Hyperslab OPerator Suite
(HOPS)](Hyperslab_OPerator_Suite_(HOPS)), and [Noesys](Noesys).

InterFormat {#InterFormat}
-----------------------------------------

[InterFormat](http://www.radio-logic.com/) is a medical image format
conversion program with both Motif and character interfaces. InterFormat
can automatically identify and convert most popular medical image
formats and write output files in many standard medical image formats,
or in formats such as netCDF that are suitable for input to leading
scientific visualization packages. InterFormat runs on UNIX
workstations; a version for OpenVMS is also available. A separate
external module for [IBM Data Explorer](#OpenDX) is available for use in
IBM Data Explorer's Visual Program Editor.
For more details about the formats handled, program features, and
pricing, see the Radio-Logic web site at
[\<http://www.radio-logic.com\>](http://www.radio-logic.com).

IRIS Explorer Module {#IRIS_Explorer_Module}
-----------------------------------------------------------

The Atmospheric and Oceanic Sciences Group at the National Center for
Supercomputing Applications (NCSA) and the Mesoscale Dynamics and
Precipitation Branch at NASA-Goddard Space Flight Center have developed
the NCSA PATHFINDER module set for [IRIS
Explorer](http://www.nag.co.uk:70/1h/Welcome_IEC). Two of the modules,
[ReadDFG](http://redrock.ncsa.uiuc.edu/PATHFINDER/pathrel2/explorer/ReadDFG/ReadDFG.html)
(to output Grids), and
[ReadDF](http://redrock.ncsa.uiuc.edu/PATHFINDER/pathrel2/explorer/ReadDF/ReadDF.html)
(to output Lattices) are capable of reading from NCSA HDF files,
MFHDF/3.3 files, and Unidata netCDF files. A user-friendly interface
provides control and information about the contents of the files.

For ReadDF, the format translation is handled transparently. Up to five
unique lattices may be generated from the file (as these files can
contain multiple data fields) using a single module. A variety of
dimensionalities and data types are supported also. Multiple variables
may be combined in a single lattice to generate vector data. All three
Explorer coordinate systems are supported.

With ReadDFG, user selected variables from the file are output in up to
five PATHFINDER grids. Each grid can consist of scalar data from one
variable or vector data from multiple variables. Coordinate information
from the file is also included in the grids. Any number of dimensions in
any of the Explorer coordinate types are supported.

For more information on the NCSA PATHFINDER project and other available
modules, visit the WWW/Mosaic PATHFINDER Home Page at
<http://redrock.ncsa.uiuc.edu/PATHFINDER/pathrel2/top/top.html> The
ReadDF module may be downloaded either via the WWW server or anonymous
ftp at redrock.ncsa.uiuc.edu in the /pub/PATHFINDER directory. For more
information please send email to: pathfinder@redrock.ncsa.uiuc.edu

See also the information on [DDI](#DDI) for another way to use netCDF
data with IRIS Explorer.

LeoNetCDF {#LeoNetCDF}
-------------------------------------

[LeoNetCDF](http://www.leokrut.com/leonetcdf.html) is a Windows
application (Windows95/NT and higher) for editing netCDF files. It can
display content of netCDF files in tree style control and permits
editing its parameters in a standard Windows interface environment.

Mathematica {#Mathematica}
-----------------------------------------

[Mathematica](http://www.wolfram.com/products/mathematica/index.html) is
a technical computing environment that provides advanced numerical and
symbolic computation and visualization. As of version 6, Mathematica
adds classic [netCDF
data](http://reference.wolfram.com/mathematica/ref/format/NetCDF.html)
to the many forms of data it can import, export, and visualize.

MATLAB {#MATLAB}
-------------------------------

[MATLAB](http://www.mathworks.com/products/matlab/) is an integrated
technical computing environment that combines numeric computation,
advanced graphics and visualization, and a high-level programming
language. Versions 7.7 and later of MATLAB have built-in support for
reading and writing netCDF data. MATLAB version 2012a includes the
netCDF 4.1.2 library with OPeNDAP client support turned on, so remote
access to netCDF and other data formats supported by OPeNDAP servers is
available.
For earlier versions, several freely-available software packages that
implement a MATLAB/netCDF interface are available:
[nctoolbox](#nctoolbox), [NetCDF Toolbox for MATLAB-5](#NC4ML5),
[MexEPS](#MexEPS), the [CSIRO MATLAB/netCDF interface](#CSIRO-MATLAB),
[NetCDF
reader](http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=15177&objectType=file),
and [fanmat](/software/netcdf/Contrib.html).



Noesys {#Neosys}
-------------------------------

[Noesys](http://www.rsinc.com/NOeSYS/index.cfm) is software for desktop
science data access and visualization. Available for both Windows and
Power Macintosh platforms, Noesys allows users to access, process,
organize and visualize large amounts of technical data.
Noesys can be used to:

-   Access and organize complex technical data
-   Export data objects to text and binary
-   View and edit large multidimensional data sets (up to 7D) in a
    spreadsheet-like environment
-   Manipulate and process data using
    [IDL®](http://www.exelisvis.com/ProductsServices/IDL.aspx), the
    Interactive Data Language, from Research Systems, Inc.
-   Interactively visualize column, matrix, and volumetric data sets
-   Image global datasets as various map projections
-   Create various projections from partial data or partial projections
    from global data (Windows only)
-   View and Edit HDF-EOS grid object data
-   Subset datasets and data tables with a GUI dialog
-   Change and save the number format of datasets and data table fields
-   Drag and Drop HDF objects between files to organize or subset files
-   Attach text annotations directly to the data file
-   Add new data objects to files and create hierarchical groups
-   Edit or create new color palettes
-   Generate publication-quality graphics for data presentation

Noesys has an interface to IDL®, allowing data to move back and forth
between Noesys and IDL with the click of a mouse. Noesys includes the
visual data analysis tools, Transform, T3D and Plot, for menu driven
plotting, rendering, and image analysis. Noesys can import HDF, HDF-EOS,
netCDF, ASCII, Binary, DTED, GeoTIFF, SDTS, TIFF, PICT, and BMP files,
create annotations, macros, images, projections and color palettes
specific to the data and save it the result as an HDF file. Noesys also
includes an HDF-EOS Grid Editor. Noesys runs on Windows 95/98 & NT and
Power Macintosh OS. More details and information about ordering Noesys
are available from
[\<http://www.rsinc.com/NOeSYS/index.cfm\>](http://www.rsinc.com/NOeSYS/index.cfm).

Origin {#Origin}
-------------------------------

Ryan Toomey reports:

Our website is <http://www.originlab.com/>

A general description of Origin: Origin includes a suite of features
that cater to the needs of scientists and engineers alike. Multi-sheet
workbooks, publication-quality graphics, and standardized analysis tools
provide a tightly integrated workspace for you to import data, create
and annotate graphs, explore and analyze data, and publish your work. To
ensure that Origin meets your data analysis requirements, intuitive
tools for advanced statistics, regression, nonlinear curve fitting,
signal processing, image processing and peak analysis are built-in.
Since any analysis operation can be set to automatically recalculate,
you can reuse your projects as templates for future work, thereby
simplifying your daily routine.

A general description of OriginPro: OriginPro offers all of the features
of Origin plus extended analysis tools for statistics, 3D fitting, image
processing and signal processing.

A general description of OriginLab Corporation: "OriginLab Corporation
produces professional data analysis and graphing software for scientists
and engineers. Our products are designed to be easy-to-use, yet have the
power and versatility to provide for the most demanding user."

PPLUS {#PPLUS}
-----------------------------

[Plot-Plus (PPLUS)](http://dwd6.home.mindspring.com/) is a general
purpose scientific graphics package, which is used in several PMEL
applications. It will read most standard ascii or binary files, as well
as netCDF file format, which used by the TOGA-TAO Project and the EPIC
system for management display and analysis. PPLUS is an interactive,
command driven, scientific graphics package which includes features such
as Mercator projection, Polar Stereographic projection, color or gray
scale area-fill contour plotting, and support for many devices:
X-windows, PostScript, HP, Tektronix, and others. This powerful and
flexible package recognizes netCDF data format, and it can extract axis
lables and graph titles from the data files. The user can customize a
plots, or combine several plots into a composite. Plots are of
publication quality. The PPLUS graphics package is used for all the TAO
workstation displays, including the animations. The animations are
created by generating a PPLUS plot for each frame, transforming the
PPLUS metacode files into HDF format with the PPLUS m2hdf filter, and
then displaying the resulting bit maps as an animation with the
XDataSlice utility, which is freely available on Internet from the
National Center for Supercomputing Applications, at
anonymous@ftp.ncsa.uiuc.edu (141.142.20.50). There is also a new m2gif
utility which produces GIF files from PPLUS metacode files.
PPLUS is supported for most Unix systems and for VAX/VMS, and is in use
at many oceanographic institutes in the US (e.g., (PMEL, Harvard, WHOI,
Scripps, NCAR, NASA, University of Rhode Island, University of Oregon,
Texas A&M...) and also internationally (Japan, Germany, Australia,
Korea...).

Plot Plus is now available at no charge. It does require licensing on a
per computer basis, but the license is at no cost. For more information
about licensing, see
[http://dwd6.home.mindspring.com/pplus_license.html/](http://dwd6.home.mindspring.com/pplus_license.html);
source and documentation are available via anonymous FTP from
<ftp://ftp.halcyon.com/pub/users/dwd/pplus1_3_2.tar.gz> and
<ftp://ftp.pmel.noaa.gov/epic/manual-dir/pplus.pdf>.

        Email:      plot_plus@halcyon.com
        Postal mail:    c/o Donald Denbo
                2138 N 186th St
                Shoreline, WA 98133
        Fax and Voice:  (206) 366-0624

PV-Wave {#PV-Wave}
---------------------------------

[PV-Wave](http://www.vni.com/products/wave/index.html) is a software
environment from [Visual Numerics](http://www.vni.com/) for solving
problems requiring the application of graphics, mathematics, numerics
and statistics to data and equations.
PV-WAVE uses a fourth generation language (4GL) that analyzes and
displays data as you enter commands. PV-WAVE includes integrated
graphics, numerics, data I/O, and data management. The latest version of
PV-Wave supports data access in numerous formats, including netCDF.

See also the information on [DDI](#DDI) for another way to use netCDF
data with PV-Wave.

Slicer Dicer {#SlicerDicer}
------------------------------------------

[Slicer Dicer](http://www.slicerdicer.com/) is a volumetric data
visualization tool, currently available for Windows and under
development for other platforms. The Slicer Dicer Web site includes a
complete list of features, an on-line user's guide, and examples of
Slicer Dicer output. Visualizations features include:
-   Perspective view of data rendered on interactively selected
    orthogonal slices, oblique slices, blocks (arbitrary rectilinear
    sub-volumes), cutouts, isosurfaces, and projected volumes (projected
    maximum, minimum, maximum absolute, or minimum absolute).
-   Optional annotations: caption, axes ticks and labels (default
    "pretty" ticks, or override to place ticks where you want them),
    color legend, data-cube outline.
-   Animation modes: slices, space, time (any parametric dimension),
    transparency, oblique slice orientation, rotation. Built-in
    animation viewer supports speed and image size controls,
    single-step, forward, backward, loop, and back-and-forth modes.
-   Select color scale from 25+ built in color tables, or import from
    palette file. Any data level or range of levels can be painted with
    an arbitrary color.
-   Any data level or range of levels can be rendered as either opaque
    or transparent.

Surfer {#Surfer}
------------------------------------------

[Surfer](http://www.goldensoftware.com/products/surfer) is a 3D
visualization, contouring, and surface modeling package that runs
under Microsoft Windows. Surfer is useful for terrain modeling,
bathymetric modeling, landscape visualization, surface analysis,
contour mapping, watershed and 3D surface mapping, gridding,
volumetrics, and more. A sophisticated interpolation engine transforms
XYZ data into publication-quality maps. Surfer imports from and
exports to a multitude of file formats, including NetCDF grids.

vGeo {#vGeo}
---------------------------

[vGeo](http://www.vrco.com/products/vgeo/vgeo.html) (Virtual Global
Explorer and Observatory) is an end-user product from
[VRCO](http://www.vrco.com/) designed to import and visualize multiple
disparate data sets, including computer simulations, observed
measurements, images, model objects, and more. vGeo is available for
IRIX, Linux and Windows platforms and supports displays ranging from
desktop monitors to multi-walled projection systems. It accepts data in
a variety of formats, including netCDF, and allows the user to specify
how multiple files and variables are mapped into a data source. 3D
graphics are built from the underlying data in real-time, and the user
has interactive control of graphics, navigation, animation, and more.

VISAGE and Decimate {#VISAGE_and_Decimate}
---------------------------------------------------------

[VISAGE](http://www.crd.ge.com/esl/cgsp/projects/visage/)
(VISualization, Animation, and Graphics Environment) is a turnkey 3D
visualization system developed at General Electric Corporate Research
and Development, (Schroeder, WJ et al, "VISAGE: An Object-Oriented
Scientific Visualization System", Proceedings of Visualization \`92
Conference). VISAGE is designed to interface with a wide variety of
data, and uses netCDF as the preferred format.

VISAGE is used at GE Corporate R & D, GE Aircraft Engine, GE Canada, GE
Power Generation, as well as ETH Zurich, Switzerland, MQS In Chieti,
Italy, and Rensselaer Polytechnic Institute in Troy, New York.

GE has another application called "Decimate" that does polygon
reduction/decimation (Schroeder,WJ et al, "Decimation of Triangle
Meshes", Proceedings of SIGGRAPH \`92). This application uses netCDF as
a preferred format. Decimate is currently licensed to Cyberware, Inc.,
makers of 3D laser digitizing hardware. Decimate is currently bundled
with the scanners, and will soon be available as a commercial product.

Voyager {#Voyager}
---------------------------------

[Makai Voyager](http://voyager.makai.com/), developed by Makai Ocean
Engineering, Inc., is 3D/4D geospatial visualization software that
enables users to import, fuse, view, and analyze large earth, ocean, and
atmosphere scientific data as it is collected or simulated in a global
geo-referenced GIS platform. The key differentiator of Makai Voyager is
its level-of-detail (LOD) technology that enables users to stream big
data rapidly over a network or the web.

Features in Makai Voyager Version 1.2 include:

-   Preprocessing LiDAR, GIS, & volumetric data from common formats into
    streamable files
-   Volume rendering for large 4D (3D + time) data, such as NetCDF
-   Analysis tools and customizable graphs
-   WMS and other streamable formats

Individual or group licenses are available for Windows (32- and 64-bit),
Linux, and Mac OS X. A full-featured 30-day trial version of Makai
Voyager is [available for download](http://voyager.makai.com).
