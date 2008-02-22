# @file SConscript
# @brief build info
#
# $Header: /nfs/slac/g/glast/ground/cvs/astro/SConscript,v 1.5 2008/02/21 22:47:29 burnett Exp $

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

progEnv.Tool('astroLib')
test_astro = progEnv.Program('test_astro', listFiles(['src/test/*.cxx']))

progEnv.Tool('registerObjects', package = 'astro', libraries = [astroLib], testApps = [test_astro], includes = listFiles(['astro/*.h']))
