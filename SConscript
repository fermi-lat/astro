# -*- python -*-
# @file SConscript
# @brief build info
#
# $Id: SConscript,v 1.77 2012/06/26 00:25:28 lsrea Exp $
# Authors: T. Burnett <tburnett@u.washington.edu>
# Version: astro-03-14-00
Import('baseEnv')
Import('listFiles')
libEnv = baseEnv.Clone()
progEnv = baseEnv.Clone()

if baseEnv['PLATFORM'] == "win32":
    libEnv.Append(CCFLAGS = "/wd4554") #  warning C4554: '<<' : check operator precedence

if 'makeStatic' in baseEnv:
    libEnv.Tool('addLinkDeps', package="astro", toBuild="static")
    astroLib = libEnv.StaticLibrary('astro', listFiles(
        ['src/*.cxx', 
         'src/wcslib/*.c', 
         'src/jplephem/*.cxx',
         'src/igrf_sub/*.cxx']))

else:
    libEnv.Tool('addLinkDeps', package="astro", toBuild="shared")
    astroLib = libEnv.SharedLibrary('astro', listFiles(
        ['src/*.cxx', 
         'src/wcslib/*.c', 
         'src/jplephem/*.cxx',
         'src/igrf_sub/*.cxx']))

progEnv.Tool('astroLib')

test_astro = progEnv.Program('test_astro', listFiles(['src/test/*.cxx']))

if 'makeStatic' in baseEnv:
    progEnv.Tool('registerTargets', 
                 package = 'astro', 
                 staticLibraryCxts = [[astroLib, libEnv]], 
                 testAppCxts = [[test_astro, progEnv]],
                 data = listFiles(['data/*']),
                 includes = listFiles(['astro/*.h']))
else:
    progEnv.Tool('registerTargets', 
                 package = 'astro', 
                 libraryCxts = [[astroLib, libEnv]], 
                 testAppCxts = [[test_astro, progEnv]],
                 data = listFiles(['data/*']),
                 includes = listFiles(['astro/*.h']))




