# @file SConscript
# @brief build info
#
# $Header: /nfs/slac/g/glast/ground/cvs/astro/SConscript,v 1.4 2008/01/31 22:50:17 burnett Exp $

import glob,os

Import('baseEnv')
Import('listFiles')
libEnv = baseEnv.Clone()
progEnv = baseEnv.Clone()

if libEnv['PLATFORM'] == "win32":
  libEnv.AppendUnique(CCFLAGS = "/D_USE_MATH_DEFINES")
  libEnv.AppendUnique(CCFLAGS = "/EHsc")
  libEnv.AppendUnique(CCFLAGS = "/D_CRT_SECURE_NO_DEPRECATE")

astroLib = libEnv.SharedLibrary('astro', listFiles(
['src/*.cxx', 
 'src/healpix/*.cxx', 
 'src/healpix/*.cc', 
 'src/wcslib/*.cxx', 
 'src/wcslib/*.c', 
 'src/jplephem/*.cxx',
 'src/igrf_sub/*.cxx']))

progEnv.Tool('registerObjects', package = 'astro', libraries = [astroLib], includes = listFiles(['astro/*.h']))


