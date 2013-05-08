sources = ['main.cc',
           'spar.cc',
           '../libqform/libqform.a',
           '../liboptarith/liboptarithxx.a',
           '/usr/local/lib/libgmp.a']

libs = ['rt', 'm']

ccflags = ['-O3', '-std=c++11', '-Wall', '-Werror']

Program(target = 'spar',
        source = sources,
        CPPPATH = '..',
        CCFLAGS = ccflags,
        LIBS = libs)
