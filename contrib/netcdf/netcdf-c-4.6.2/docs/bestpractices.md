Writing NetCDF Files: Best Practices {#BestPractices}
====================================

[TOC]

Best Practices {#bp_Best_Practices}
=====================================

## Conventions {#bp_Conventions}

While netCDF is intended for "self-documenting data", it is often
necessary for data writers and readers to agree upon attribute
conventions and representations for discipline-specific data structures.
These agreements are written up as human readable documents called
***netCDF conventions***.

Use an existing Convention if possible. See the list of [registered
conventions](/software/netcdf/conventions.html).

The CF Conventions are recommended where applicable, especially for
gridded (model) datasets.

Document the convention you are using by adding the global attribute
"Conventions" to each netCDF file, for example:

This document refers to conventions for the netCDF *classic* data model.
For recommendations about conventions for the netCDF-4 *enhanced* data
model, see [Developing Conventions for
NetCDF-4](/netcdf/papers/nc4_conventions.html).

## Coordinate Systems {#bp_Coordinate-Systems}

A ***coordinate variable*** is a one-dimensional variable with the same
name as a dimension, which names the coordinate values of the dimension.
It must not have any missing data (for example, no `_FillValue` or
`missing_value` attributes) and must be strictly monotonic (values
increasing or decreasing). A two-dimensional variable of type char is a
***string-valued coordinate variable*** if it has the same name as its
first dimension, e.g.: **char time( time, time\_len);** all of its
strings must be unique. A variable's ***coordinate system*** is the set
of coordinate variables used by the variable. Coordinates that refer to
physical space are called ***spatial coordinates***, ones that refer to
physical time are called ***time coordinates***, ones that refer to
either physical space or time are called ***spatio-temporal
coordinates.***

-   Make coordinate variables for every dimension possible (except for
    string length dimensions).
-   Give each coordinate variable at least `unit` and `long_name`
    attributes to document its meaning.
-   Use an existing netCDF [Convention](#Conventions) for your
    coordinate variables, especially to identify
    spatio-temporal coordinates.
-   Use shared dimensions to indicate that two variables use the same
    coordinates along that dimension. If two variables' dimensions are
    not related, create separate dimensions for them, even if they
    happen to have the same length.

## Variable Grouping {#bp_Variable-Grouping}

You may structure the data in a netCDF file in different ways, for
example putting related parameters into a single variable by adding an
extra dimension. Standard visualization and analysis software may have
trouble breaking that data out, however. On the other extreme, it is
possible to create different variables e.g. for different vertical
levels of the same parameter. However, standard visualization and
analysis software may have trouble grouping that data back together.
Here are some guidelines for deciding how to group your data into
variables:

-   All of the data in a variable must be of the same type and should
    have the same units of measurement.
-   A variable's attributes should be applicable to all its data.
-   If possible, all of the coordinate variables should be
    spatio-temporal, with no extra dimensions.
-   Use 4D spatio-temporal coordinate systems in preference to 3D. Use
    3D spatio-temporal coordinate systems in preference to 2D.
-   Vector valued (e.g. wind) parameters are legitimate uses of extra
    dimensions. There are trade-offs between putting vectors in the same
    variables vs. putting each component of a vector in a
    different variable. Check that any visualization software you plan
    to use can deal with the structure you choose.
-   Think in terms of complete coordinate systems (especially
    spatio-temporal), and organize your data into variables accordingly.
    Variables with the same coordinate system implicitly form a group.

## Variable Attributes {#bp_Variable-Attributes}


-   For each variable where it makes sense, add a **units** attribute,
    using the [udunits](/software/udunits/index.html) conventions,
    if possible.
-   For each variable where it makes sense, add a **long\_name ****
    attribute, which is a human-readable descriptive name for
    the variable. This could be used for labeling plots, for example.

## Strings and Variables of type char {#bp_Strings-and-Variables-of-type-char}

NetCDF-3 does not have a primitive **String** type, but does have arrays
of type **char**, which are 8 bits in size. The main difference is that
Strings are variable length arrays of chars, while char arrays are fixed
length. Software written in C usually depends on Strings being zero
terminated, while software in Fortran and Java do not. Both C
(*nc\_get\_vara\_text*) and Java (*ArrayChar.getString*) libraries have
convenience routines that read char arrays and convert to Strings.

-   Do not use char type variables for numeric data, use byte type
    variables instead.
-   Consider using a global Attribute instead of a Variable to store a
    String applicable to the whole dataset.
-   When you want to store arrays of Strings, use a multidimensional
    char array. All of the Strings will be the same length.
-   There are 3 strategies for writing variable length Strings and
    zero-byte termination:
    1.  *Fortran convention*: pad with blanks and never terminate with a
        zero byte.
    2.  *C convention*: pad with zeros and always terminate with a
        zero byte.
    3.  *Java convention*: You don't need to store a trailing zero byte,
        but pad trailing unused characters with zero bytes.
-   When reading, trim zeros and blanks from the end of the char array
    and if in C, add a zero byte terminator.

## Calendar Date/Time {#bp_Calendar-Date-Time}

Time as a fundamental unit means a time interval, measured in seconds. A
Calendar date/time is a specific instance in real, physical time. Dates
are specified as an interval from some ***reference time*** e.g. "days
elapsed since Greenwich mean noon on 1 January 4713 BCE". The reference
time implies a system of counting time called a ***calendar*** (e.g.
Gregorian calendar) and a textual representation (e.g. [ISO
8601](http://www.cl.cam.ac.uk/%7Emgk25/iso-time.html)).

There are two strategies for storing a date/time into a netCDF variable.
One is to encode it as a numeric value and a unit that includes the
reference time, e.g. "seconds since 2001-1-1 0:0:0" or"days since
2001-1-1 0:0:0" . The other is to store it as a String using a standard
encoding and Calendar. The former is more compact if you have more than
one date, and makes it easier to compute intervals between two dates.

Unidata's [udunits](/software/udunits/) package provides a convenient
way to implement the first strategy. It uses the ISO 8601 encoding and a
hybrid Gregorian/Julian calendar, but udunits does not support use of
other Calendars or encodings for the reference time. However the ncdump
"-T" option can display numeric times that use udunits (and optionally
climate calendars) as ISO 8601 strings that are easy for humans to
interpret.

-   If your data uses real, physical time that is well represented using
    the Gregorian/Julian calendar, encode it as an interval from a
    reference time, and add a units attribute which uses a
    udunits-compatible time unit. If the data assumes one of the
    non-standard calendars mentioned in the CF Conventions, specify that
    with a Calendar attribute. Readers can then use the udunits package
    to manipulate or format the date values, and the ncdump utility can
    display them with either numeric or string representation.
-   If your data uses a calendar not supported by the CF Conventions,
    make it compatible with existing date manipulation packages if
    possible (for example, java.text.SimpleDateFormat).
-   Add multiple sets of time encodings if necessary to allow different
    readers to work as well as possible.\

## Unsigned Data {#bp_Unsigned-Data}

NetCDF-3 does not have unsigned integer primitive types.

-   To be completely safe with unknown readers, widen the data type, or
    use floating point.
-   You can use the corresponding signed types to store unsigned data
    only if all client programs know how to interpret this correctly.
-   A new proposed convention is to create a variable attribute
    `_Unsigned = "true"` to indicate that integer data should be treated
    as unsigned.

## Packed Data Values {#bp_Packed-Data-Values}

Packed data is stored in a netCDF file by limiting precision and using a
smaller data type than the original data, for example, packing
double-precision (64-bit) values into short (16-bit) integers. The
C-based netCDF libraries do not do the packing and unpacking. (The
[netCDF Java library](/software/netcdf-java/) will do automatic
unpacking when the
[VariableEnhanced](/software/netcdf-java/v4.1/javadocAll/ucar/nc2/dataset/VariableEnhanced.html)
Interface is used. For details see
[EnhancedScaleMissing](/software/netcdf-java/v4.1/javadocAll/ucar/nc2/dataset/EnhanceScaleMissing.html)).

-   Each variable with packed data has two attributes called
    **scale\_factor** and **add\_offset**, so that the packed data may
    be read and unpacked using the formula:

    > ***unpacked\_data\_value = packed\_data\_value \* scale\_factor +
    > add\_offset***

-   The type of the stored variable is the packed data type, typically
    byte, short or int.
-   The type of the scale\_factor and add\_offset attributes should be
    the type that you want the unpacked data to be, typically float
    or double.
-   To avoid introducing a bias into the unpacked values due to
    truncation when packing, the data provider should round to the
    nearest integer rather than just truncating towards zero before
    writing the data:

    > ***packed\_data\_value = nint((unpacked\_data\_value -
    > add\_offset) / scale\_factor)***

Depending on whether the packed data values are intended to be
interpreted by the reader as signed or unsigned integers, there are
alternative ways for the data provider to compute the *scale\_factor*
and *add\_offset* attributes. In either case, the formulas above apply
for unpacking and packing the data.

A conventional way to indicate whether a byte, short, or int variable is
meant to be interpreted as unsigned, even for the netCDF-3 classic model
that has no external unsigned integer type, is by providing the special
variable attribute `_Unsigned` with value `"true"`. However, most
existing data for which packed values are intended to be interpreted as
unsigned are stored without this attribute, so readers must be aware of
packing assumptions in this case. In the enhanced netCDF-4 data model,
packed integers may be declared to be of the appropriate unsigned type.

Let *n* be the number of bits in the packed type, and assume *dataMin*
and *dataMax* are the minimum and maximum values that will be used for a
variable to be packed.

-   If the packed values are intended to be interpreted as signed
    integers (the default assumption for classic model data), you may
    use:

    > *scale\_factor =(dataMax - dataMin) / (2^n^ - 1)*

    > *add\_offset = dataMin + 2^n\\ -\\ 1^ \* scale\_factor*

-   If the packed values are intended to be interpreted as unsigned (for
    example, when read in the C interface using the `nc_get_var_uchar()`
    function), use:

    > *scale\_factor =(dataMax - dataMin) / (2^n^ - 1)*

    > *add\_offset = dataMin*

-   In either the signed or unsigned case, an alternate formula may be
    used for the add\_offset and scale\_factor packing parameters that
    reserves a packed value for a special value, such as an indicator of
    missing data. For example, to reserve the minimum packed value
    (-2^n\\ -\\ 1^) for use as a special value in the case of signed
    packed values:

    > *scale\_factor =(dataMax - dataMin) / (2^n^ - 2)*

    > *add\_offset = (dataMax + dataMin) / 2*

-   If the packed values are unsigned, then the analogous formula that
    reserves 0 as the packed form of a special value would be:

    > *scale\_factor =(dataMax - dataMin) / (2^n^ - 2)*

    > *add\_offset = dataMin - scale\_factor*

-   Example, packing 32-bit floats into 16-bit shorts:

            variables:
              short data( z, y, x);
                data:scale_offset = 34.02f;
                data:add_offset = 1.54f;

-   The `units` attribute applies to unpacked values.

## Missing Data Values {#bp_Missing-Data-Values}

***Missing data*** is a general name for data values that are invalid,
never written, or missing. The netCDF library itself does not handle
these values in any special way, except that the value of a `_FillValue`
attribute, if any, is used in pre-filling unwritten data. (The
Java-netCDF library will assist in recognizing these values when
reading, see class **VariableStandardized**).

-   Default fill values for each type are available in the C-based
    interfaces, and are defined in the appropriate header files. For
    example, in the C interface, NC\_FILL\_FLOAT and NC\_FILL\_DOUBLE
    are numbers near 9.9692e+36 that are returned when you try to read
    values that were never written. Writing, reading, and testing for
    equality with these default fill values works portably on the
    platforms on which netCDF has been tested.
-   The `_FillValue` attribute should have the same data type as the
    variable it describes. If the variable is packed using
    `scale_factor` and `add_offset` attributes, the `_FillValue`
    attribute should have the data type of the packed data.
-   Another way of indicating missing values for real type data is to
    store an IEEE **NaN** floating point value. The advantage of this is
    that any computation using a NaN results in a NaN. Client software
    must know to look for NaNs, however, and detection of NaNs is
    tricky, since any comparison with a NaN is required to return
    *false*.
    -   In Java, you can use **Double.NaN** and **Float.NaN** constants.
    -   In many C compilers, you can generate a NaN value using **double
        nan = 0.0 / 0.0;**
-   Alternatively or in addition, set the **valid\_range** attribute for
    each variable that uses missing values, and make sure all valid data
    is within that range, and all missing or invalid data is outside of
    that range. Again, the client software must recognize and make use
    of this information. Example:

            variables:
              float data( z, y, x);
                data:valid_range = -999.0f, 999.0f;


    If the variable is packed using `scale_factor` and `add_offset`
    attributes, the `valid_range` attribute should have the data type of
    the packed data.

    If the variable is unsigned the `valid_range` values should be
    widened if needed and stored as unsigned integers.

## Miscellaneous tips {#bp_Miscellaneous-tips}

-   To define a file whose structure is known in advance, write a CDL
    file and create the netCDF file using
    [ncgen](/cgi-bin/man-cgi?ncgen). Then write the data into the netCDF
    file using your program. This is typically much easier than
    programming all of the create calls yourself.
-   For the netCDF classic or 64-bit-offset formats, it's possible to
    reserve extra space in the file when it is created so that you can
    later add additional attributes or non-record variables without
    copying all the data. (This is not necessary for netCDF-4 files,
    because metadata can be added efficiently to such files.) See the [C
    man-page reference documentation](/cgi-bin/man-cgi?netcdf+-s3) (or
    the [Fortran reference documentation](/cgi-bin/man-cgi?netcdf+-s3f))
    for `nc__create` and `nc__enddef` (`nf__create` and `nf__enddef`
    for Fortran) for more details on reserving extra space in
    the header.

## Spelling netCDF: Best Practices {#bp_Spelling-netCDF-Best-Practices}

There are only 3 correct spellings of "netCDF":

1.  **netCDF:** The original spelling of the name of the data model,
    API, and format. The acronym stands for network Common Data Form
    (not Format), and the "CDF" part was capitalized in part to pay
    homage to the NASA "CDF" data model which the netCDF data
    model extended.
2.  **netcdf:** Used in certain file names, such as:

             #include <netcdf.h>  

3.  **NetCDF**: Used in titles and at the beginning of sentences, where
    "netCDF" is awkward or violates style guidelines.

All other forms, and most especially "Netcdf", are considered vulgar and
a sign of ill-breeding or misspent youth, analogous to the egregious but
common misspelling "JAVA" used by those who are new to the language or
who mistakenly think it is an acronym.\
