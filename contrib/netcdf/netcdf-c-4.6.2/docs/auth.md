NetCDF Authorization Support
======================================

<!-- double header is needed to workaround doxygen bug -->

NetCDF Authorization Support {#auth}
====================================

[TOC]

## Introduction {#auth_intro}

netCDF can support user authorization using the facilities provided by the curl
library. This includes basic password authentication as well as
certificate-based authorization.

At the moment, this document only applies to DAP2 and DAP4 access
because they are (for now) the only parts of the netCDF-C library
that uses libcurl.

With some exceptions (e.g. see the section on <a href="#REDIR">redirection</a>)
The libcurl authorization mechanisms can be accessed in two ways

1. Inserting the username and password into the url, or
2. Accessing information from a so-called _rc_ file named either
   `.daprc` or `.dodsrc`

## URL-Based Authentication {#auth_url}

For simple password based authentication, it is possible to
directly insert the username and the password into a url in this form.

    http://username:password@host/...

This username and password will be used if the server asks for
authentication. Note that only simple password authentication
is supported in this format.

Specifically note that [redirection-based](#REDIR)
authorization may work with this but it is a security risk.
This is because the username and password
may be sent to each server in the redirection chain.

Note also that the `user:password` form may contain characters that must be
escaped. See the <a href="#USERPWDESCAPE">password escaping</a> section to see
how to properly escape the user and password.

## RC File Authentication {#auth_dodsrc}
The netcdf library supports an _rc_ file mechanism to allow the passing
of a number of parameters to libnetcdf and libcurl.
Locating the _rc_ file is a multi-step process.

### Search Order

The file must be called one of the following names:
".daprc" or ".dodsrc".
If both ".daprc" and ".dodsrc" exist, then
the ".daprc" file will take precedence.

It is strongly suggested that you pick one of the two names
and use it always. Otherwise you may observe unexpected results
when the netcdf-c library finds one that you did not intend.

The search for an _rc_ file looks in the following places in this order.

1. Check for the environment variable named _DAPRCFILE_.
   This will specify the full path for the _rc_ file
   (not just the containing directory).
2. Search the current working directory (`./`) looking
   for (in order) .daprc or .dodsrc.
3. Search the HOME directory (`$HOME`) looking
   for (in order) .daprc or .dodsrc. The HOME environment
   variable is used to define the directory in which to search.

It is strongly suggested that you pick a uniform location
and use it always. Otherwise you may observe unexpected results
when the netcdf-c library get an rc file you did not expect.

### RC File Format

The rc file format is a series of lines of the general form:

    [<host:port>]<key>=<value>

where the bracket-enclosed host:port is optional.

### URL Constrained RC File Entries

Each line of the rc file can begin with
a host+port enclosed in square brackets.
The form is "host:port".
If the port is not specified
then the form is just "host".
The reason that more of the url is not used is that
libcurl's authorization grain is not any finer than host level.

Examples.

    [remotetest.unidata.ucar.edu]HTTP.VERBOSE=1

or

    [fake.ucar.edu:9090]HTTP.VERBOSE=0

If the url request from, say, the _netcdf_open_ method
has a host+port matching one of the prefixes in the rc file, then
the corresponding entry will be used, otherwise ignored.
This means that an entry with a matching host+port will take
precedence over an entry without a host+port.

For example, the URL

    http://remotetest.unidata.ucar.edu/thredds/dodsC/testdata/testData.nc

will have HTTP.VERBOSE set to 1 because its host matches the example above.

Similarly,

    http://fake.ucar.edu:9090/dts/test.01

will have HTTP.VERBOSE set to 0 because its host+port matches the example above.

## Authorization-Related Keys {#auth_keys}

The currently defined set of authorization-related keys are as follows.
The second column is the affected curl_easy_setopt option(s), if any
(see reference #1).
<table>
<tr><th>Key</th><th>Affected curl_easy_setopt Options</th><th>Notes</th>
<tr><td>HTTP.COOKIEJAR</td><td>CURLOPT_COOKIEJAR</td>
<tr><td>HTTP.COOKIEFILE</td><td>CURLOPT_COOKIEJAR</td><td>Alias for CURLOPT_COOKIEJAR</td>
<tr><td>HTTP.PROXY.SERVER</td><td>CURLOPT_PROXY, CURLOPT_PROXYPORT, CURLOPT_PROXYUSERPWD</td>
<tr><td>HTTP.PROXY_SERVER</td><td>CURLOPT_PROXY, CURLOPT_PROXYPORT, CURLOPT_PROXYUSERPWD</td><td>Decprecated: use HTTP.PROXY.SERVER</td>
<tr><td>HTTP.SSL.CERTIFICATE</td><td>CURLOPT_SSLCERT</td>
<tr><td>HTTP.SSL.KEY</td><td>CURLOPT_SSLKEY</td>
<tr><td>HTTP.SSL.KEYPASSWORD</td><td>CURLOPT_KEYPASSWORD</td>
<tr><td>HTTP.SSL.CAINFO</td><td>CURLOPT_CAINFO</td>
<tr><td>HTTP.SSL.CAPATH</td><td>CURLOPT_CAPATH</td>
<tr><td>HTTP.SSL.VERIFYPEER</td><td>CURLOPT_SSL_VERIFYPEER</td>
<tr><td>HTTP.SSL.VALIDATE</td><td>CURLOPT_SSL_VERIFYPEER, CURLOPT_SSL_VERIFYHOST</td>
<tr><td>HTTP.CREDENTIALS.USERPASSWORD</td><td>CURLOPT_USERPASSWORD</td>
<tr><td>HTTP.CREDENTIALS.USERNAME</td><td>CURLOPT_USERNAME</td>
<tr><td>HTTP.CREDENTIALS.PASSWORD</td><td>CURLOPT_PASSWORD</td>
<tr><td>HTTP.NETRC</td><td>N.A.</td><td>Specify path of the .netrc file</td>
</table>

### Password Authentication

The key
HTTP.CREDENTIALS.USERPASSWORD
can be used to set the simple password authentication.
This is an alternative to setting it in the url.
The value must be of the form "username:password".
See the <a href="#USERPWDESCAPE">password escaping</a> section
to see how this value must escape certain characters.
Also see <a href="#REDIR">redirection authorization</a>
for important additional information.

The pair of keys
HTTP.CREDENTIALS.USERNAME and HTTP.CREDENTIALS.PASSWORD
can be used as an alternative to HTTP.CREDENTIALS.USERPASSWORD
to set the simple password authentication.
If present, they take precedence over HTTP.CREDENTIALS.USERPASSWORD.
The values do not need to be escaped.
See <a href="#REDIR">redirection authorization</a>
for important additional information.

### Cookie Jar

The HTTP.COOKIEJAR key
specifies the name of file from which
to read cookies (CURLOPT_COOKIEJAR) and also
the file into which to store cookies (CURLOPT_COOKIEFILE).
The same value is used for both CURLOPT values.
It defaults to in-memory storage.
See [redirection authorization](#REDIR)
for important additional information.

### Certificate Authentication

HTTP.SSL.CERTIFICATE
specifies a file path for a file containing a PEM cerficate.
This is typically used for client-side authentication.

HTTP.SSL.KEY is essentially the same as HTTP.SSL.CERTIFICATE
and should always have the same value.

HTTP.SSL.KEYPASSWORD
specifies the password for accessing the HTTP.SSL.CERTIFICAT/HTTP.SSL.key file.

HTTP.SSL.CAPATH
specifies the path to a directory containing
trusted certificates for validating server certificates.
See reference #2 for more info.

HTTP.SSL.VALIDATE
is a boolean (1/0) value that if true (1)
specifies that the client should verify the server's presented certificate.

HTTP.PROXY.SERVER
specifies the url for accessing the proxy:
e.g. *http://[username:password@]host[:port]*

HTTP.PROXY_SERVER
deprecated; use HTTP.PROXY.SERVER

HTTP.NETRC
specifies the absolute path of the .netrc file.
See [redirection authorization](#REDIR)
for information about using .netrc.

## Password Escaping {#auth_userpwdescape}

With current password rules, it is is not unlikely that the password
will contain characters that need to be escaped. Similarly, the user
may contain characters such as '@' that need to be escaped. To support this,
it is assumed that all occurrences of `user:password` use URL (i.e. %%XX)
escaping for at least the characters in the table below.

The minimum set of characters that must be escaped depends on the location.
If the user+pwd is embedded in the URL, then '@' and ':' __must__ be escaped.
If the user+pwd is the value for 
the HTTP.CREDENTIALS.USERPASSWORD key in the _rc_ file, then
':' __must__ be escaped.
Escaping should __not__ be used in the `.netrc` file nor in
HTTP.CREDENTIALS.USERNAME or HTTPCREDENTIALS.PASSWORD.

The relevant escape codes are as follows.
<table>
<tr><th>Character</th><th>Escaped Form</th>
<tr><td>'@'</td><td>%40</td>
<tr><td>':'</td><td>%3a</td>
</table>
Additional characters can be escaped if desired.

## Redirection-Based Authentication {#auth_redir}

Some sites provide authentication by using a third party site
to do the authentication. Examples include ESG, URS, RDA, and most oauth2-based
systems.

The process is usually as follows.

1. The client contacts the server of interest (SOI), the actual data provider
using, typically _http_ protocol.
2. The SOI sends a redirect to the client to connect to the e.g. URS system
using the _https_ protocol (note the use of _https_ instead of _http_).
3. The client authenticates with URS.
4. URS sends a redirect (with authorization information) to send
the client back to the SOI to actually obtain the data.

It turns out that libcurl, by default, uses the password in the
`.daprc` file (or from the url) for all connections that request
a password.  This causes problems because only the the specific
redirected connection is the one that actually requires the password.
This is where the `.netrc` file comes in. Libcurl will use `.netrc`
for the redirected connection. It is possible to cause libcurl
to use the `.daprc` password always, but this introduces a
security hole because it may send the initial user+pwd to every
server in the redirection chain.
In summary, if you are using redirection, then you are
''strongly'' encouraged to create a `.netrc` file to hold the
password for the site to which the redirection is sent.

The format of this `.netrc` file will contain lines that
typically look like this.

    machine mmmmmm login xxxxxx password yyyyyy

where the machine, mmmmmm, is the hostname of the machine to
which the client is redirected for authorization, and the
login and password are those needed to authenticate on that machine.

The location of the `.netrc` file can be specified by
putting the following line in your `.daprc`/`.dodsrc` file.

    HTTP.NETRC=<path to netrc file>

If not specified, then libcurl will look first in the current
directory, and then in the HOME directory.

One final note. In using this, you MUST
to specify a real file in the file system to act as the
cookie jar file (HTTP.COOKIEJAR) so that the
redirect site can properly pass back authorization information.

## Client-Side Certificates {#auth_clientcerts}

Some systems, notably ESG (Earth System Grid), requires
the use of client-side certificates, as well as being
[re-direction based](#REDIR).
This requires setting the following entries:

- HTTP.COOKIEJAR &mdash; a file path for storing cookies across re-direction.
- HTTP.NETRC &mdash; the path to the netrc file.
- HTTP.SSL.CERTIFICATE &mdash; the file path for the client side certificate file.
- HTTP.SSL.KEY &mdash; this should have the same value as HTTP.SSL.CERTIFICATE.
- HTTP.SSL.CAPATH &mdash; the path to a "certificates" directory.
- HTTP.SSL.VALIDATE &mdash; force validation of the server certificate.

Note that the first two are there to support re-direction based authentication.

## References

1. https://curl.haxx.se/libcurl/c/curl_easy_setopt.html
2. https://curl.haxx.se/docs/ssl-compared.html

## Appendix A. All RC-File Keys {#auth_allkeys}

For completeness, this is the list of all rc-file keys.
If this documentation is out of date with respect to the actual code,
the code is definitive.
<table>
<tr><th>Key</th><th>curl_easy_setopt Option</th>
<tr valign="top"><td>HTTP.DEFLATE</td><td>CUROPT_DEFLATE<br>with value "deflate,gzip"</td>
<tr><td>HTTP.VERBOSE</td><td>CUROPT_VERBOSE</td>
<tr><td>HTTP.TIMEOUT</td><td>CUROPT_TIMEOUT</td>
<tr><td>HTTP.USERAGENT</td><td>CUROPT_USERAGENT</td>
<tr><td>HTTP.COOKIEJAR</td><td>CUROPT_COOKIEJAR</td>
<tr><td>HTTP.COOKIE_JAR</td><td>CUROPT_COOKIEJAR</td>
<tr valign="top"><td>HTTP.PROXY.SERVER</td><td>CURLOPT_PROXY,<br>CURLOPT_PROXYPORT,<br>CURLOPT_PROXYUSERPWD</td>
<tr valign="top"><td>HTTP.PROXY_SERVER</td><td>CURLOPT_PROXY,<br>CURLOPT_PROXYPORT,<br>CURLOPT_PROXYUSERPWD</td>
<tr><td>HTTP.SSL.CERTIFICATE</td><td>CUROPT_SSLCERT</td>
<tr><td>HTTP.SSL.KEY</td><td>CUROPT_SSLKEY</td>
<tr><td>HTTP.SSL.KEYPASSWORD</td><td>CUROPT_KEYPASSWORD</td>
<tr><td>HTTP.SSL.CAINFO</td><td>CUROPT_CAINFO</td>
<tr><td>HTTP.SSL.CAPATH</td><td>CUROPT_CAPATH</td>
<tr><td>HTTP.SSL.VERIFYPEER</td><td>CUROPT_SSL_VERIFYPEER</td>
<tr><td>HTTP.CREDENTIALS.USERPASSWORD</td><td>CUROPT_USERPASSWORD</td>
<tr><td>HTTP.CREDENTIALS.USERNAME</td><td>CUROPT_USERNAME</td>
<tr><td>HTTP.CREDENTIALS.PASSWORD</td><td>CUROPT_PASSWORD</td>
<tr><td>HTTP.NETRC</td><td>CURLOPT_NETRC,CURLOPT_NETRC_FILE</td>
</table>

## Appendix B. URS Access in Detail {#auth_ursdetail}

It is possible to use the NASA Earthdata Login System (URS)
with netcdf by using using the process specified in the
[redirection based authorization section](#REDIR).
In order to access URS controlled datasets, however, it is necessary to
register as a user with NASA at this website (subject to change):

    https://uat.urs.earthdata.nasa.gov/

## Appendix C. ESG Access in Detail {#auth_esgdetail}

It is possible to access Earth Systems Grid (ESG) datasets
from ESG servers through the netCDF API using the techniques
described in the section on [Client-Side Certificates](#CLIENTCERTS).

In order to access ESG datasets, however, it is necessary to
register as a user with ESG and to setup your environment
so that proper authentication is established between an netcdf
client program and the ESG data server.  Specifically, it
is necessary to use what is called "client-side keys" to
enable this authentication. Normally, when a client accesses
a server in a secure fashion (using "https"), the server
provides an authentication certificate to the client.
With client-side keys, the client must also provide a
certificate to the server so that the server can know with
whom it is communicating. Note that this section is subject
to change as ESG changes its procedures.

The netcdf library uses the _curl_ library and it is that
underlying library that must be properly configured.

### Terminology

The key elements for client-side keys requires the constructions of
two "stores" on the client side.

* Keystore - a repository to hold the client side key.
* Truststore - a repository to hold a chain of certificates
that can be used to validate the certificate
sent by the server to the client.

The server actually has a similar set of stores, but the client
need not be concerned with those.

### Initial Steps

The first step is to obtain authorization from ESG.
Note that this information may evolve over time, and
may be out of date.
This discussion is in terms of BADC and NCSA. You will need
to substitute as necessary.

1. Register at http://badc.nerc.ac.uk/register
   to obtain access to badc and to obtain an openid,
   which will looks something like:
   <pre>https://ceda.ac.uk/openid/Firstname.Lastname</pre>

2. Ask BADC for access to whatever datasets are of interest.

3. Obtain short term credentials at
   _http://grid.ncsa.illinois.edu/myproxy/MyProxyLogon/_
   You will need to download and run the MyProxyLogon program.
   This will create a keyfile in, typically, the directory ".globus".
   The keyfile will have a name similar to this: "x509up_u13615"
   The other elements in ".globus" are certificates to use in
   validating the certificate your client gets from the server.

4. Obtain the program source ImportKey.java
   from this location: _http://www.agentbob.info/agentbob/79-AB.html_
   (read the whole page, it will help you understand the remaining steps).

### Building the KeyStore

You will have to modify the keyfile in the previous step
and then create a keystore and install the key and a certificate.
The commands are these:

    openssl pkcs8 -topk8 -nocrypt -in x509up_u13615 -inform PEM -out key.der -outform DER
    openssl x509 -in x509up_u13615 -inform PEM -out cert.der -outform DER
    java -classpath <path to ImportKey.class> -Dkeypassword="<password>" -Dkeystore=./<keystorefilename> key.der cert.der

Note, the file names "key.der" and "cert.der" can be whatever you choose.
It is probably best to leave the .der extension, though.

### Building the TrustStore

Building the truststore is a bit tricky because as provided, the
certificates in ".globus" need some massaging. See the script below
for the details. The primary command is this, which is executed for every
certificate, c, in globus. It sticks the certificate into the file
named "truststore"

    keytool -trustcacerts -storepass "password" -v -keystore "truststore"  -importcert -file "${c}"

### Running the C Client

Refer to the section on [Client-Side Certificates](#CLIENTCERTS).
The keys specified there  must be set in the rc file to support ESG access.

- HTTP.COOKIEJAR=~/.dods_cookies
- HTTP.NETRC=~/.netrc
- HTTP.SSL.CERTIFICATE=~/esgkeystore
- HTTP.SSL.KEY=~/esgkeystore
- HTTP.SSL.CAPATH=~/.globus
- HTTP.SSL.VALIDATE=1

Of course, the file paths above are suggestions only;
you can modify as needed.
The HTTP.SSL.CERTIFICATE and HTTP.SSL.KEY
entries should have same value, which is the file path for the
certificate produced by MyProxyLogon.  The HTTP.SSL.CAPATH entry
should be the path to the "certificates" directory produced by
MyProxyLogon.

As noted, ESG also uses re-direction based authentication.
So, when it receives an initial connection from a client, it
redirects to a separate authentication server. When that
server has authenticated the client, it redirects back to
the original url to complete the request.

### Script for creating Stores

The following script shows in detail how to actually construct the key
and trust stores. It is specific to the format of the globus file
as it was when ESG support was first added. It may have changed
since then, in which case, you will need to seek some help
in fixing this script. It would help if you communicated
what you changed to the author so this document can be updated.

    #!/bin/sh -x
    KEYSTORE="esgkeystore"
    TRUSTSTORE="esgtruststore"
    GLOBUS="globus"
    TRUSTROOT="certificates"
    CERT="x509up_u13615"
    TRUSTROOTPATH="$GLOBUS/$TRUSTROOT"
    CERTFILE="$GLOBUS/$CERT"
    PWD="password"

    D="-Dglobus=$GLOBUS"
    CCP="bcprov-jdk16-145.jar"
    CP="./build:${CCP}"
    JAR="myproxy.jar"

    # Initialize needed directories
    rm -fr build
    mkdir build
    rm -fr $GLOBUS
    mkdir $GLOBUS
    rm -f $KEYSTORE
    rm -f $TRUSTSTORE

    # Compile MyProxyCmd and ImportKey
    javac -d ./build -classpath "$CCP" *.java
    javac -d ./build ImportKey.java

    # Execute MyProxyCmd
    java -cp "$CP myproxy.MyProxyCmd

    # Build the keystore
    openssl pkcs8 -topk8 -nocrypt -in $CERTFILE -inform PEM -out key.der -outform DER
    openssl x509 -in $CERTFILE -inform PEM -out cert.der -outform DER
    java -Dkeypassword=$PWD -Dkeystore=./${KEYSTORE} -cp ./build ImportKey key.der cert.der

    # Clean up the certificates in the globus directory
    for c in ${TRUSTROOTPATH}/*.0 ; do
        alias=`basename $c .0`
        sed -e '0,/---/d' <$c >/tmp/${alias}
        echo "-----BEGIN CERTIFICATE-----" >$c
        cat /tmp/${alias} >>$c
    done

    # Build the truststore
    for c in ${TRUSTROOTPATH}/*.0 ; do
        alias=`basename $c .0`
        echo "adding: $TRUSTROOTPATH/${c}"
        echo "alias: $alias"
        yes | keytool -trustcacerts -storepass "$PWD" -v -keystore ./$TRUSTSTORE -alias $alias -importcert -file "${c}"
    done
    exit

## Point of Contact

__Author__: Dennis Heimbigner<br>
__Email__: dmh at ucar dot edu
__Initial Version__: 11/21/2014<br>
__Last Revised__: 08/24/2017

