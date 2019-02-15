NetCDF Credits {#credits}
==============

[Unidata](http://www.unidata.ucar.edu/) is sponsored by the [National
Science Foundation](http://www.nsf.gov/) and managed by the [University
Corporation for Atmospheric Research](http://www2.ucar.edu/).

The NASA CDF data model, to which netCDF owes much, was developed by
Michael Gough and Lloyd Treinish. Joe Fahle designed a C version for a
CDF-like interface and discussions with Joe provided much of the
inspiration for the original netCDF C interface. The netCDF C library
and first version of the Java interface were written by Glenn Davis. The
nctest test suite, ncdump, ncgen, the original C++ interface, and most
of the original User's Guides were written by Russ Rew. The Fortran77
test code was written by Cathy Cormack. The configure-based installation
system, vars and varm array access implementation, Fortran77 interface,
and original perl interface were written by Steve Emmerson. The netCDF-3
interface design and the nc\_test exhaustive test code were developed by
Glenn Davis and Harvey Davies. The Fortran77 interface uses Burkhard
Burow's cfortran.h package. John Caron wrote the subsequent netCDF Java
implementations. The Fortran90 interface and User's Guide were developed
by Robert Pincus. Ed Hartnett updated and simplified the configure-based
installation, enhanced Windows support, refactored the documentation
while converting it into texinfo, and is the primary developer for
netCDF-4. Dennis Heimbigner wrote the netCDF-4 version of ncgen, the C
OPeNDAP client, the dispatch layer, and the implementation of diskless
files. The nccopy utility was added by Russ Rew. Lynton Appel developed
the C++ implementation for netCDF-4. Ward Fisher overhauled netCDF release-engineering, developed a new build-and-test framework using CMake, virtualization, and container technologies, moved sources to GitHub, developed a Windows/Microsoft Visual Studio port, refactored documentation for improved web access, and merged all the documentation into the sources for maintaining with Markdown and Doxygen.

The following people have contributed related software, bug reports,
fixes, valuable suggestions, and other kinds of useful support:

Jennifer Adams, Bob Albrecht, Chuck Alexander, Dave Allured, Ethan
Alpert, Chris Anderson, Ayal Anis, Harald Anlauf, Lynton Appel,
Sylwester Arabas, Stephen Armstrong, Ben Auer, Phil Austin, Fedor Baart,
Eric Bachalo, Jason Bacon, Larry Baker, Sandy Ballard, Matthew Banta,
Christopher Bartz, Sourish Basu, Mike Berkley, Stephen Bespalko, Ingo
Bethke, Sherman Beus, Sachin Kumar Bhate, John A. Biddiscombe, Lorenzo
Bigagli, Mark Borges, Nicola Botta, Kenneth P. Bowman, Bill Boyd, Mark
Bradford, Bernward Bretthauer, Paul A. Bristow, Roy Britten, Dave Brown,
Alexander Bruhns, Ryan Cabell, Peter Cao, Jed O. Caplan, Glenn Carver,
Tom Cavin, Morrell Chance, Susan C. Cherniss, Jason E. Christy, Gerardo
Cisneros, Alain Coat, Carlie J. Coats, Jr., Antonio S. Cofiño, Tony
Conrad, Jon Corbet, Alexandru Corlan, Gus Correa, Jim Cowie, Alex
Crosby, Arlindo da Silva, Chris Dallimore, Rick Danielson, Alan Dawes,
Donald W. Denbo, Charles R. Denham, Arnaud Desitter, Sue Dettling, Steve
Diggs, Martin Dix, Michael Dixon, Jacek Dlugosz, Alastair Doherty,
Charles Doutriaux, Bob Drach, Ludovic Drouineau, Huaiyu Duan, Patrice
Dumas, Paul J. Durack, Frank Dzaak, Brian Eaton, Harry Edmon, Veit
Eitner, Lee Elson, Mario Emmenlauer, Takeshi Enomoto, Ata Etemadi,
Constantinos Evangelinos, John Evans, Joe Fahle, Gabor Fichtinger, Glenn
Flierl, Connor J. Flynn, Shanna-Shaye Forbes, Anne Fouilloux,
Jean-Francois Foccroulle, Mike Folk, David Forrest, David W. Forslund,
Ben Foster, Masaki Fukuda, Dave Fulker, Ryo Furue, James Gallagher, Jose
Luis Garcia, Rao Garimella, Bear Giles, Tom Glaess, Peter Gleckler,
Christoph Gohlke, Paul Goodman, Nath Gopalaswamy, André Gosselin, Udo
Grabowski, Gary Granger, Jonathan Gregory, Markus Gross, Lionel Guez,
Patrick Guio, Mark Hadfield, Magnus Hagdorn, Paul Hamer, Steve Hankin,
Petr Hanousek, Christopher Harrop, Bill Hart, David Hassell, Ros
Hatcher, Mark Hedley, Kate Hedstrom, Charles Hemphill, Barron Henderson,
Joerg Henrichs, Olaf Heudecker, Thijs Heus, Donn Hines, Konrad Hinsen,
Yuan Ho, Kari Hoijarvi, Leigh Holcombe, Tim Holt, Toshinobu Hondo,
Takeshi Horinouchi, Chris Houck, Wei Huang, Matt Huddleston, Matt
Hughes, Tim Hume, Doug Hunt, Nathanael Hübbe, L. F. Hwang, Alan Imerito,
Jouk Jansen, Steve Jardin, Rimvydas Jasinskas, Harry Jenter, Susan
Jesuroga, Patrick Jöckel, Tomas Johannesson, Peter Gylling Jørgensen,
Junchang Ju. Narita Kazumi, Maxwell Kelley, John Kemp, Jamie
Kettleborough, Nikolay Khabarov, Constantine Khroulev, Steve Kirby,
Heiko Klein, Luis Kornblueh, Frank Kreienkamp, Jeff Kuehn, Julian
Kunkel, Jarle Ladstein, V. Lakshmanan, Bruce Langdon, Éric Larouche,
Stephen Leak, Benoit Lecocq, Tom LeFebvre, JF Le Fillátre, Angel Li,
Jianwei Li, Zhi Liang, Wei-keng Liao, Si Liu, Rick Light, Brian Lincoln,
Keith Lindsay, Fei Liu, Josep Llodrà, Jeffery W. Long, Dave Lucas,
Valerio Luccio, Lifeng Luo, John Lillibridge, Steve Luzmoor, Lawrence
Lyjak, Rich Lysakowski, Sergey Malyshev, Len Makin, Harry Mangalam,
Ansley Manke, David Mansbach, Jim Mansbridge, Andreas Manschke, Chris
Marquardt, Marinna Martini, William C. Mattison, Craig Mattocks, Sean
McBride, Mike McCarrick, Bill McKie, Ron Melton, Roy Mendelssohn, Pavel
Michna, Barb Mihalas, Henry LeRoy Miller Jr., Philip Miller, Rakesh
Mithal, Masahiro Miiyaki, Kengo Miyamoto, Tushar Mohan, Christine C.
Molling, Skip Montanaro, Thomas L. Moore, Paidemwoyo Munhutu, Stefano
Nativi, Gottfried Necker, Peter Neelin, Erik Noble, Michael Nolta, Bill
Noon, Enda O'Brien, Thomas Orgis, Dave Osburn, Dan Packman, Mark Payne,
Simon Paech, Doug Palmer, Gabor Papp, Stephen Pascoe, Morten Pedersen,
Louise Perkins, Michael D Perryman, Hartmut Peters, Ron Pfaff, David
Pierce, Alexander Pletzer, Philippe Poilbarbe, Dierk Polzin, Orion
Poplawski, Jacob Weismann Poulsen, Ken Prada, Greg Rappa, Dave Raymond,
Michael Redetzky, Rene Redler, Mark Reyes, Doug Reynolds, Jose Agustín
García Reynoso, Mike Rilee, Mark Rivers, Randolph Roesler, Mike Romberg,
Mathis Rosenhauer, Cédric Roux, Suzanne T. Rupert, Stephen Sachs,
Toshihiro Sakakima, Eric Salathe, Sean Patrick Santos, Matthew H.
Savoie, Marie Schall, Remko Scharroo, Brian Schlining, Nico Schlömer.
Dan Schmitt, Robert B. Schmunk, Larry A. Schoof, Rich Schramm, William
J. Schroeder, Karen Schuchardt, Uwe Schulzweida, Keith Searight, Johann
Schmitz, Andreas Schwab, Guntram Seiß, John Sheldon, Sergei Shibaev,
Masato Shiotani, Michael Shopsin, Richard P. Signell, Steve Simpson, Joe
Sirott, Greg Sjaardema, Dirk Slawinski, Cathy Smith, Neil R. Smith,
Ruben Smits, Peter P. Smolka, Nancy Soreide, Hudson Souza, Gunter
Spranz, Richard Stallman, Bjorn Stevens, Reto Stöckli, Bob Swanson, John
Tanski, Karl Taylor, Jason Thaxter, Kevin W. Thomas, Olivier Titaud,
Jonathan Tomshine, Mark Tracy, Philippe Tulkens, Warren Turkal, Tom
Umeda, Joe VanAndel, Paul van Delst, Gerald van der Grijn, Luuk van
Dijk, Richard van Hees, Martin van Driel, Jànos Vègh, Jos Verdoold,
Matthieu Verstraete, Pedro Vicente, Lykle Voort, Hailin Yan, Lianqing
Yu, Bernhard Wagner, Thomas Wainwright, Stephen Walker, Ya-Qiang Wang,
John C. Warner, Chris Webster, Richard Weed, Paul Wessel, Jeffrey S.
Whitaker, Carsten Wieczorrek, Gerry Wiener, Ralf Wildenhues, David
Wilensky, Hartmut Wilhelms, Gareth Williams, Florian Wobbe, David
Wojtowicz, Jeff Wong, Randy Zagar, Charlie Zender, Remik Ziemlinski.

Development and implementation of netCDF is supported by the National
Science Foundation, Unidata's primary sponsor. Development of the
netCDF-4 interface was initially funded by NASA's Earth Science
Technology Office. Addition of OPeNDAP client support to the netCDF
library is based on work supported by the National Science Foundation.
Any opinions, findings and conclusions or recomendations expressed in
this material are those of the authors and do not necessarily reflect
the views of the sponsoring organizations.
