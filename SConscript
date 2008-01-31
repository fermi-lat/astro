# @file SConscript
# @brief build info
#
# $Header: /nfs/slac/g/glast/ground/cvs/astro/SConscript,v 1.3 2008/01/02 22:29:26 burnett Exp $

import glob,os

Import('baseEnv')
Import('listFiles')
libEnv = baseEnv.Clone()
progEnv = baseEnv.Clone()

astroLib = libEnv.SharedLibrary('astro', listFiles(
['src/*.cxx', 
 'src/healpix/*.cxx', 
 'src/healpix/*.cc', 
 'src/wcslib/*.cxx', 
 'src/wcslib/*.c', 
 'src/jplephem/*.cxx',
 'src/igrf_sub/*.cxx']))

progEnv.Tool('registerObjects', package = 'astro', libraries = [astroLib], includes = listFiles(['astro/*.h']))
