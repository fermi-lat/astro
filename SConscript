# -*- python -*-
# @file SConscript
# @brief build info
#
# $Id: SConscript,v 1.30 2009/02/06 21:30:04 glastrm Exp $
# Authors: T. Burnett <tburnett@u.washington.edu>
# Version: astro-03-06-01
Import('baseEnv')
Import('listFiles')
libEnv = baseEnv.Clone()
progEnv = baseEnv.Clone()

if baseEnv['PLATFORM'] == "win32":
    libEnv.Append(CCFLAGS = "/wd4554") #  warning C4554: '<<' : check operator precedence

libEnv.Tool('astroLib', depsOnly = 1)
astroLib = libEnv.SharedLibrary('astro', listFiles(
['src/*.cxx', 
 'src/wcslib/*.c', 
 'src/jplephem/*.cxx',
 'src/igrf_sub/*.cxx']))

progEnv.Tool('astroLib')
test_astro = progEnv.Program('test_astro', listFiles(['src/test/*.cxx']))

progEnv.Tool('registerObjects', 
             package = 'astro', 
	     libraries = [astroLib], 
	     testApps = [test_astro], 
	     includes = listFiles(['astro/*.h']))




