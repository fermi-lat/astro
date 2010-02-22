# -*- python -*-
# @file SConscript
# @brief build info
#
# $Id: SConscript,v 1.46 2010/02/18 00:49:53 jrb Exp $
# Authors: T. Burnett <tburnett@u.washington.edu>
# Version: astro-03-09-00
Import('baseEnv')
Import('listFiles')
libEnv = baseEnv.Clone()
progEnv = baseEnv.Clone()

if baseEnv['PLATFORM'] == "win32":
    libEnv.Append(CCFLAGS = "/wd4554") #  warning C4554: '<<' : check operator precedence

if baseEnv['PLATFORM'] == "win32":
    libEnv.Tool('astroLib', depsOnly = 1)
    
astroLib = libEnv.SharedLibrary('astro', listFiles(
['src/*.cxx', 
 'src/wcslib/*.c', 
 'src/jplephem/*.cxx',
 'src/igrf_sub/*.cxx']))

progEnv.Tool('astroLib')
test_astro = progEnv.Program('test_astro', listFiles(['src/test/*.cxx']))

progEnv.Tool('registerTargets', 
             package = 'astro', 
	     libraryCxts = [[astroLib, libEnv]], 
	     testAppCxts = [[test_astro, progEnv]],
             data = listFiles(['data/*']),
	     includes = listFiles(['astro/*.h']))




