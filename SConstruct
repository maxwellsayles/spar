cflags = ['-O3', '-Wall', '-Werror']
cxxflags = ['-O3', '-std=c++11', '-Wall', '-Werror']

library_sources = ['spar.cc',
                   'primorials/s64_theory_opt.c',
                   'primorials/s128_theory_opt.c']

StaticLibrary(target='spar',
              source=library_sources,
              CPPPATH = '..',
              CFLAGS = cflags,
              CXXFLAGS = cxxflags)
              
sources = ['main.cc',
           'libspar.a',
           '../libqform/libqform.a',
           '../liboptarith/liboptarithxx.a',
           '/usr/local/lib/libgmp.a']

libs = ['rt', 'm']

Program(target = 'spar',
        source = sources,
        CPPPATH = '..',
        CFLAGS = cflags,
        CXXFLAGS = cxxflags,	
        LIBS = libs)
