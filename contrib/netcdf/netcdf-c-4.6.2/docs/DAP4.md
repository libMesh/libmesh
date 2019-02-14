# DAP4 Introduction {#dap4_intro}

Beginning with netCDF version 4.5.0, optional support is provided for
accessing data through servers supporting the DAP4 protocol.

DAP4 support is enabled if the _--enable-dap__ option
is used with _./configure_. If DAP4 support is enabled, then
a usable version of _libcurl_ must be specified
using the _LDFLAGS_ environment variable (similar to the way
that the _HDF5_ libraries are referenced).
Refer to the installation manual for details.
By default DAP4 support is enabled if _libcurl_ is found.
DAP4 support can be disabled using the _--disable-dap_.

DAP4 uses a data model that is by design similar to the netCDF enhanced
(netCDF-4) data model.
Generically, the DAP4 meta-data is encoded textually in a _DMR_
object. For detailed information about DAP4, refer to the
DAP4 specification
https://docs.opendap.org/index.php/OPULS_Development#DAP4_Specification

# Accessing DAP4 Data {#dap4_accessing_data}

In order to access a DAP4 data source through the netCDF API, the
file name normally used is replaced with a URL with a specific
format. The URL is composed of three parts.
-  URL - this is a standard form URL such as
   http://remotetest.unidata.ucar.edu/d4ts/test.01

-  Constraints - these are suffixed to the URL in the query part
   of the url (“?\<constraint\>”. The structure of the constraint
   is somewhat complicated; and DAP4 specification 
   should be consulted. The interaction of DAP4
   constraints with netCDF is complex and at the moment requires an
   understanding of how DAP4 is translated to netCDF.

- Client parameters - these may be specified in either of
  two ways.  The older, deprecated form prefixes text to the
  front of the url and is of the the general form [\<name>]
  or [\<name>=value].  Examples include [show=fetch] and
  [noprefetch].  The newer, preferred form prefixes the
  parameters to the end of the url using the semi-standard '#'
  format: e.g. http://....#show=fetch&noprefetch.

It is possible to see what the translation does to a particular
DAP4 data source by examining the DMR source through a web
browser and then examining the translation using the _ncdump -h_
command to see the netCDF Classic translation.

For example, if a web browser is given the following, this URL
will return the DMR in XML format for the specified dataset.
````
http://149.165.169.123:8080/d4ts/testfiles/test_one_var.nc.dmr.xml
````

Then by using the following ncdump command, it is possible to see the
equivalent netCDF enhanced translation.
````
ncdump -h http://149.165.169.123:8080/d4ts/testfiles/test_one_var.nc#dap4
````
Note the use of the "dap4" fragment added at the end. This tells the netCDF
library to use the DAP4 protocol instead of the default DAP2 protocol.

The DMR output from the web server should look like this.
````
<Dataset name="test_one_var.nc" dapVersion="4.0" dmrVersion="1.0">
    <Int32 name="t"/>
    <Attribute name="_DAP4_Little_Endian" type="UInt8">
        <Value value="1"/>
    </Attribute>
</Dataset>
````

The output from ncdump should look like this.
````
netcdf test_one_var {
variables:
	int t ;

// global attributes:
		:_DAP4_Little_Endian = 1UB ;
}
````

# DAP4 to NetCDF Translation Rules {#dap4_to_netcdf}

The netCDF library DAP4 code translates the DAP4 data model
into the netCDF enhanced (netCDF-4) data model.

## netCDF-4 Translation Rules {#dap4_nc42_trans_rules}

For illustrative purposes, the following example DMR will be used.
````
<Dataset name="test_groups1.nc" dapVersion="4.0" dmrVersion="1.0">
    <Dimension name="dim1" size="5"/>
    <Attribute name="_DAP4_Little_Endian" type="UInt8"><Value value="1"/></Attribute>
    <Group name="g">
        <Dimension name="dim2" size="3"/>
        <Group name="h">
            <Dimension name="dim3" size="7"/>
            <Int32 name="v1"><Dim name="/dim1"/></Int32>
            <Float32 name="v2"><Dim name="/g/dim2"/></Float32>
        </Group>
        <Group name="i">
            <Dimension name="dim3" size="7"/>
            <Int32 name="v1"><Dim name="/dim1"/></Int32>
            <Float32 name="v3"><Dim name="/g/i/dim3"/></Float32>
        </Group>
    </Group>
</Dataset>
````
This will translate (via ncdump) into this.
````
netcdf test_groups1 {
dimensions:
	dim1 = 5 ;

// global attributes:
		:_DAP4_Little_Endian = 1UB ;

group: g {
  dimensions:
  	dim2 = 3 ;

  group: h {
    dimensions:
    	dim3 = 7 ;
    variables:
    	int v1(dim1) ;
    	float v2(dim2) ;
    } // group h

  group: i {
    dimensions:
    	dim3 = 7 ;
    variables:
    	int v1(dim1) ;
    	float v3(dim3) ;
    } // group i
  } // group g
}
````

## Variable Definition {#dap4_var2_def}

The set of netCDF variables is derived from the fields with primitive
base types as they occur in Sequences, Grids, and Structures. The
field names are modified to be fully qualified initially. For the
above, the set of variables are as follows. The coordinate variables
within grids are left out in order to mimic the behavior of libnc-dap4.
````
    f1
    S1.f11
    S1.FS2.f1
    S1.FS2.f2
    S2.G1.temp
    S2.G2.G2
    lat
    lon
````

## DAP4 Reserved Keywords {#dap4_reserved_keywords}

In the OPeNDAP DAP4 protocol, there are a number of reserved keywords.  These keywords are case insensitive and if you use one as a netCDF variable name, you may encounter odd behavior such as case changes (depending on the client DDS/DAS parser).  The list of reserved keywords as used by the netCDF-C library parser are as follows:

- alias
- array
- attributes
- byte
- dataset
- error
- float32
- float64
- grid
- int16
- int32
- maps
- sequence
- string
- structure
- uint16
- uint32
- url
- code
- message
- program_type
- program


## Variable Dimension Translation {#dap4_var_dim_trans}

A variable's rank is determined from three sources.
- The variable has the dimensions associated with the field it
represents (e.g. S1.FS2.f1[3] in the above example).
- The variable inherits the dimensions associated with any containing
structure that has a rank greater than zero. These dimensions precede
those of case 1. Thus, we have in our example, f1[2][3], where the
first dimension comes from the containing Structure FS2[2].
- The variable's set of dimensions are altered if any of its
containers is a DAP4 DDS Sequence. This is discussed more fully below.

If the type of the netCDF variable is char, then an extra string
dimension is added as the last dimension.

## Dimension translation {#dap4_dim2_trans}

For dimensions, the rules are as follows.

Fields in dimensioned structures inherit the dimension of the
structure; thus the above list would have the following dimensioned
variables.
````
        S1.FS2.f1 -> S1.FS2.f1[2][3]
        S1.FS2.f2 -> S1.FS2.f2[2]
        S2.G1.temp -> S2.G1.temp[lat=2][lon=2]
        S2.G1.lat -> S2.G1.lat[lat=2]
        S2.G1.lon -> S2.G1.lon[lon=2]
        S2.G2.G2 -> S2.G2.lon[lat=2][lon=2]
        S2.G2.lat -> S2.G2.lat[lat=2]
        S2.G2.lon -> S2.G2.lon[lon=2]
        lat -> lat[lat=2]
        lon -> lon[lon=2]
````

Collect all of the dimension specifications from the DDS, both named
and anonymous (unnamed) For each unique anonymous dimension with value
NN create a netCDF dimension of the form "XX_\<i\>=NN", where XX is the
fully qualified name of the variable and i is the i'th (inherited)
dimension of the array where the anonymous dimension occurs. For our
example, this would create the following dimensions.
````
        S1.FS2.f1_0 = 2 ;
        S1.FS2.f1_1 = 3 ;
        S1.FS2.f2_0 = 2 ;
        S2.G2.lat_0 = 2 ;
        S2.G2.lon_0 = 2 ;
````

If however, the anonymous dimension is the single dimension of a MAP
vector in a Grid then the dimension is given the same name as the map
vector This leads to the following.
````
        S2.G2.lat_0 -> S2.G2.lat
        S2.G2.lon_0 -> S2.G2.lon
````

For each unique named dimension "<name>=NN", create a netCDF dimension
of the form "<name>=NN", where name has the qualifications removed. If
this leads to duplicates (i.e. same name and same value), then the
duplicates are ignored. This produces the following.
````
        S2.G2.lat -> lat
        S2.G2.lon -> lon
````

Note that this produces duplicates that will be ignored later.

At this point the only dimensions left to process should be named
dimensions with the same name as some dimension from step number 3,
but with a different value. For those dimensions create a dimension of
the form "<name>M=NN" where M is a counter starting at 1. The example
has no instances of this.

Finally and if needed, define a single UNLIMITED dimension named
"unlimited" with value zero. Unlimited will be used to handle certain
kinds of DAP4 sequences (see below).

This leads to the following set of dimensions.
````
dimensions:
  unlimited = UNLIMITED;
  lat = 2 ;
  lon = 2 ;
  S1.FS2.f1_0 = 2 ;
  S1.FS2.f1_1 = 3 ;
  S1.FS2.f2_0 = 2 ;
````

## Variable Name Translation {#dap4_var_name_trans}

The steps for variable name translation are as follows.

Take the set of variables captured above. Thus for the above DDS, the
following fields would be collected.
````
        f1
        S1.f11
        S1.FS2.f1
        S1.FS2.f2
        S2.G1.temp
        S2.G2.G2
        lat
        lon
````

All grid array variables are renamed to be the same as the containing
grid and the grid prefix is removed. In the above DDS, this results in
the following changes.
````
        G1.temp -> G1
        G2.G2 -> G2
````

It is important to note that this process could produce duplicate
variables (i.e. with the same name); in that case they are all assumed
to have the same content and the duplicates are ignored. If it turns
out that the duplicates have different content, then the translation
will not detect this. YOU HAVE BEEN WARNED.

The final netCDF-3 schema (minus attributes) is then as follows.
````
netcdf t {
dimensions:
        unlimited = UNLIMITED ;
        lat = 2 ;
        lon = 2 ;
        S1.FS2.f1_0 = 2 ;
        S1.FS2.f1_1 = 3 ;
        S1.FS2.f2_0 = 2 ;
variables:
        int f1 ;
        int lat(lat) ;
        int lon(lon) ;
        int S1.f11 ;
	int S1.FS2.f1(S1.FS2.f1_0, S1.FS2.f1_1) ;
        int S1.FS2.f2(S1_FS2_f2_0) ;
        float S2.G1(lat, lon) ;
        float G2(lat, lon) ;
}
````

In practice, the unlimited dimension is dropped because it is unused.

There are differences with the original libnc-dap4 here because
libnc-dap4 technically was incorrect. The original would have said
this, for example.
````
int S1.FS2.f1(lat, lat) ;
````

Note that this is incorrect because it dimensions S1.FS2.f1(2,2)
rather than S1.FS2.f1(2,3).

## Translating DAP4 DDS Sequences {#dap4_translation}

Any variable (as determined above) that is contained directly or
indirectly by a Sequence is subject to revision of its rank using the
following rules.

Let the variable be contained in Sequence Q1, where Q1 is the
innermost containing sequence. If Q1 is itself contained (directly or
indirectly) in a sequence, or Q1 is contained (again directly or
indirectly) in a structure that has rank greater than 0, then the
variable will have an initial UNLIMITED dimension. Further, all
dimensions coming from "above" and including (in the containment
sense) the innermost Sequence, Q1, will be removed and replaced by
that single UNLIMITED dimension. The size associated with that
UNLIMITED is zero, which means that its contents are inaccessible
through the netCDF-3 API. Again, this differs from libnc-dap4, which
leaves out such variables. Again, however, this difference is backward
compatible.

If the variable is contained in a single Sequence (i.e. not nested)
and all containing structures have rank 0, then the variable will have
an initial dimension whose size is the record count for that
Sequence. The name of the new dimension will be the name of the
Sequence.

Consider this example.
````
Dataset {
  Structure {
    Sequence {
      Int32 f1[3];
      Int32 f2;
    } SQ1;
  } S1[2];
  Sequence {
    Structure {
      Int32 x1[7];
    } S2[5];
  } Q2;
} D;
````

The corresponding netCDF-3 translation is pretty much as follows (the
value for dimension Q2 may differ).
````
dimensions:
    unlimited = UNLIMITED ; // (0 currently)
    S1.SQ1.f1_0 = 2 ;
    S1.SQ1.f1_1 = 3 ;
    S1.SQ1.f2_0 = 2 ;
    Q2.S2.x1_0 = 5 ;
    Q2.S2.x1_1 = 7 ;
    Q2 = 5 ;
variables:
    int S1.SQ1.f1(unlimited, S1.SQ1.f1_1) ;
    int S1.SQ1.f2(unlimited) ;
    int Q2.S2.x1(Q2, Q2.S2.x1_0, Q2.S2.x1_1) ;
````

Note that for example S1.SQ1.f1_0 is not actually used because it has
been folded into the unlimited dimension.

Note that for sequences without a leading unlimited dimension, there
is a performance cost because the translation code has to walk the
data to determine how many records are associated with the
sequence. Since libnc-dap4 did essentially the same thing, it can be
assumed that the cost is not prohibitive.

# Caching {#dap4_dap4_caching}

In an effort to provide better performance for some access patterns,
client-side caching of data is available. The default is no caching,
but it may be enabled by prefixing the URL with the parameter "cache".

Caching operates basically as follows.

When a URL is first accessed using _nc_open()_, netCDF automatically
does a pre-fetch of selected variables. These include all variables
smaller than a specified (and user definable) size. This allows, for
example, quick access to coordinate variables. This can be suppressed
with the parameter "noprefetch".

Whenever a request is made using some variant of the _nc_get_var()_ API
procedures, the complete variable is fetched and stored in the cache
as a new cache entry. Subsequence requests for any part of that
variable will access the cache entry to obtain the data.

The cache may become too full, either because there are too many
entries or because it is taking up too much disk space. In this case
cache entries are purged until the cache size limits are reached. The
cache purge algorithm is LRU (least recently used) so that variables
that are repeatedly referenced will tend to stay in the cache.

The cache is completely purged when _nc_close()_ is invoked.

In order to decide if you should enable caching, you will need to have
some understanding of the access patterns of your program.

The ncdump program always dumps one or more whole variables so it
turns on caching.

If your program accesses only parts of a number of variables, then
caching should probably not be used since fetching whole variables
will probably slow down your program for no purpose.

Unfortunately, caching is currently an all or nothing proposition, so
for more complex access patterns, the decision to cache or not may not
have an obvious answer. Probably a good rule of thumb is to avoid
caching initially and later turn it on to see its effect on
performance.

# Defined Client Parameters {#dap4_dap4_defined_params}

Currently, a limited set of client parameters is
recognized. Parameters not listed here are ignored, but no error is
signalled. All names are case insensitive.

Parameter Name Legal Values Semantics
- "log" | "log=<file>" - Turn on logging and send the log output to
  the specified file. If no file is specified, then log output is sent
  to standard error.
- "show=... das|dds|url" - This causes information to appear as
  specific global attributes. The currently recognized tags are "dds"
  to display the underlying DDS, "das" similarly, and "url" to display
  the url used to retrieve the data. This parameter may be specified
  multiple times (e.g. “show=dds&show=url”).
- "show=fetch" - This parameter causes the netCDF code to log a copy
  of the complete url for every HTTP get request. If logging is
  enabled, then this can be helpful in checking to see the access
  behavior of the netCDF code.
- "stringlength=NN" - Specify the default string length to use for
  string dimensions. The default is 64. The name "maxstrlen" is an
  alias for "stringlength".
- "stringlength_\<var\>=NN" - Specify the default string length to use
  for a string dimension for the specified variable. The default is
  64. The name "maxstrlen_\<var\>" is an alias for "stringlength_\<var\>".
- "cache" - This enables caching.
- "nocache" - This disbles caching.
- "cachelimit=NN" - Specify the maximum amount of space allowed for
  the cache.
- "cachecount=NN" - Specify the maximum number of entries in the
  cache.
- "prefetch" - This enables prefetch of small variables (default).
- "noprefetch" - This disables prefetch of small variables.
- "fillmismatch" - This enables _FillValue/Variable type mismatch.
- "nofillmismatch" - This disables _FillValue/Variable type mismatch (default).

# Notes on Debugging OPeNDAP Access {#dap4_dap4_debug}

The OPeNDAP support makes use of the logging facility of the
underlying oc system (see http://www.OPeNDAP.org/oc).
Note that this is currently separate from the
existing netCDF logging facility. Turning on this logging can
sometimes give important information. Logging can be enabled by
using the client parameter "log" or "log=filename",
where the first case will send log output to standard error and the
second will send log output to the specified file.

Users should also be aware that if one is
accessing data over an NFS mount, one may see some .nfsxxxxx files;
those can be ignored.

## HTTP Configuration. {#dap4_http2_config}

Limited support for configuring the http connection is provided via
parameters in the “.dodsrc” configuration file. The relevant .dodsrc file is
located by first looking in the current working directory, and if not
found, then looking in the directory specified by the “$HOME”
environment variable.

Entries in the .dodsrc file are of the form:
````
     ['['<url>']']<key>=<value>
````

That is, it consists of a key name and value pair and optionally
preceded by a url enclosed in square brackets.

For given KEY and URL strings, the value chosen is as follows:

If URL is null, then look for the .dodsrc entry that has no url prefix
and whose key is same as the KEY for which we are looking.

If the URL is not null, then look for all the .dodsrc entries that
have a url, URL1, say, and for which URL1 has the same host and port
as URL. All parts of the url's except host and port are ignored.
For example, if URL = http//x.y/a, then it will match
entries of the form
_[http//x.y/a]KEY=VALUE_ or _[http//x.y/b]KEY=VALUE_.
It will not match an entry of the form _[http//x.y:8080]KEY=VALUE
because the second has a port number (8080) different than the URL.
Finally from the set so constructed, choose the first matching entry.

Currently, the supported set of keys (with descriptions) are as
follows.

1. HTTP.VERBOSE  
        Type: boolean ("1"/"0")  
        Description: Produce verbose output, especially using SSL.  
        Related CURL Flags: CURLOPT_VERBOSE  
1. HTTP.DEFLATE  
        Type: boolean ("1"/"0")  
        Description: Allow use of compression by the server.  
        Related CURL Flags: CURLOPT_ENCODING  
1. HTTP.COOKIEJAR  
        Type: String representing file path  
        Description: Specify the name of file into which to store cookies. Defaults to in-memory storage.  
        Related CURL Flags:CURLOPT_COOKIEJAR  
1. HTTP.CREDENTIALS.USER  
        Type: String representing user name  
        Description: Specify the user name for Digest and Basic authentication.  
        Related CURL Flags:  
1. HTTP.CREDENTIALS.PASSWORD  
        Type: String representing password  
        Type: boolean ("1"/"0")  
        Description: Specify the password for Digest and Basic authentication.  
        Related CURL Flags:  
1. HTTP.SSL.CERTIFICATE  
        Type: String representing file path  
        Description: Path to a file containing a PEM cerficate.  
        Related CURL Flags: CURLOPT_CERT  
1. HTTP.SSL.KEY  
        Type: String representing file path  
        Description: Same as HTTP.SSL.CERTIFICATE, and should usually have the same value.  
        Related CURL Flags: CURLOPT_SSLKEY  
1. HTTP.SSL.KEYPASSWORD  
        Type: String representing password  
        Description: Password for accessing the HTTP.SSL.KEY/HTTP.SSL.CERTIFICATE  
        Related CURL Flags: CURLOPT_KEYPASSWORD  
1. HTTP.SSL.CAPATH  
        Type: String representing directory  
        Description: Path to a directory containing trusted certificates for validating server certificates.  
        Related CURL Flags: CURLOPT_CAPATH  
1. HTTP.SSL.VALIDATE  
        Type: boolean ("1"/"0")  
        Description: Cause the client to verify the server's presented certificate.  
        Related CURL Flags: CURLOPT_SSL_VERIFYPEER, CURLOPT_SSL_VERIFYHOST  
1. HTTP.TIMEOUT  
        Type: String ("dddddd")  
        Description: Specify the maximum time in seconds that you allow the http transfer operation to take.  
        Related CURL Flags: CURLOPT_TIMEOUT, CURLOPT_NOSIGNAL  
1. HTTP.PROXY_SERVER  
        Type: String representing url to access the proxy: (e.g.http://[username:password@]host[:port])  
        Description: Specify the needed information for accessing a proxy.  
        Related CURL Flags: CURLOPT_PROXY, CURLOPT_PROXYHOST, CURLOPT_PROXYUSERPWD  
1. HTTP.READ.BUFFERSIZE  
        Type: String ("dddddd")  
        Description: Specify the the internal buffer size for curl reads.  
        Related CURL Flags: CURLOPT_BUFFERSIZE, CURL_MAX_WRITE_SIZE (16kB),  
                            CURL_MAX_READ_SIZE (512kB).  
  
1. HTTP.KEEPALIVE  
        Type: String ("on|n/m")  
        Description: Specify that TCP KEEPALIVE should be enabled and that the associated idle wait time is n and that the associated repeat interval is m. If the value is of the form is the string "on", then turn on keepalive, but do not set idle or interval.  
        Related CURL Flags: CURLOPT_TCP_KEEPALIVE, CURLOPT_TCP_KEEPIDLE,  
                            CURLOPT_TCP_KEEPINTVL.  

The related curl flags line indicates the curl flags modified by this
key. See the libcurl documentation of the _curl_easy_setopt()_ function
for more detail (http://curl.haxx.se/libcurl/c/curl_easy_setopt.html).

For ESG client side key support, the following entries must be specified:
````
HTTP.SSL.VALIDATE
HTTP.COOKIEJAR
HTTP.SSL.CERTIFICATE
HTTP.SSL.KEY
HTTP.SSL.CAPATH
````

Additionally, for ESG, the _HTTP.SSL.CERTIFICATE_ and _HTTP.SSL.KEY_
entries should have same value, which is the file path for the
certificate produced by MyProxyLogon. The HTTP.SSL.CAPATH entry should
be the path to the "certificates" directory produced by MyProxyLogon.

# Point of Contact {#dap4_poc}

__Author__: Dennis Heimbigner<br>
__Email__: dmh at ucar dot edu<br>
__Initial Version__: 6/5/2017<br>
__Last Revised__: 9/25/2018
