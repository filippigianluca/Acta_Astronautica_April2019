 #    Copyright (c) 2010-2013, Delft University of Technology
 #    All rights reserved.
 #
 #    Redistribution and use in source and binary forms, with or without modification, are
 #    permitted provided that the following conditions are met:
 #      - Redistributions of source code must retain the above copyright notice, this list of
 #        conditions and the following disclaimer.
 #      - Redistributions in binary form must reproduce the above copyright notice, this list of
 #        conditions and the following disclaimer in the documentation and/or other materials
 #        provided with the distribution.
 #      - Neither the name of the Delft University of Technology nor the names of its contributors
 #        may be used to endorse or promote products derived from this software without specific
 #        prior written permission.
 #
 #    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 #    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 #    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 #    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 #    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 #    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 #    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 #    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 #    OF THE POSSIBILITY OF SUCH DAMAGE.
 #
 #    Changelog
 #      YYMMDD    Author            Comment
 #      151102          Niccolo Gastaldello    File for shaping functions of Spherical Shaping for GTOC2
 #
 #    References
 #
 #    Notes

# Set the source files.
set(SHAPING_SOURCES
  "${SRCROOT}${SHAPINGDIR}/basefuncderiv.cpp"
  "${SRCROOT}${SHAPINGDIR}/sphericalshaping.cpp"
  "${SRCROOT}${SHAPINGDIR}/timefunc.cpp"
  "${SRCROOT}${SHAPINGDIR}/approximatePlanets.cpp"
  "${SRCROOT}${SHAPINGDIR}/getEphemeris.cpp"
  "${SRCROOT}${SHAPINGDIR}/m2theta.cpp"
  "${SRCROOT}${SHAPINGDIR}/addexcessvelocity.cpp"
)

# Set the header files.
set(SHAPING_HEADERS
  "${SRCROOT}${SHAPINGDIR}/basefuncderiv.h"
  "${SRCROOT}${SHAPINGDIR}/sphericalshaping.h"
  "${SRCROOT}${SHAPINGDIR}/timefunc.h"
  "${SRCROOT}${SHAPINGDIR}/loaddata.h"
  "${SRCROOT}${SHAPINGDIR}/approximatePlanets.h"
  "${SRCROOT}${SHAPINGDIR}/approximatePlanetsEcliptic.h"
  "${SRCROOT}${SHAPINGDIR}/ephemerisBase.h"
  "${SRCROOT}${SHAPINGDIR}/getEphemeris.h"
  "${SRCROOT}${SHAPINGDIR}/m2theta.h"
  "${SRCROOT}${SHAPINGDIR}/addexcessvelocity.h"
)

# Add static libraries, second line only if to be used later on outside this application.
add_library(functions_shaping STATIC ${SHAPING_SOURCES} ${SHAPING_HEADERS})
