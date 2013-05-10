sources = ['main.cc',
	   'spar.cc',
	   'primorials/s64_theory_opt.c',
	   'primorials/s128_theory_opt.c',
           '../libqform/libqform.a',
           '../liboptarith/liboptarithxx.a',
           '/usr/local/lib/libgmp.a']

libs = ['rt', 'm']

cflags = ['-O3', '-Wall', '-Werror']
cxxflags = ['-O3', '-std=c++11', '-Wall', '-Werror']

Program(target = 'spar',
        source = sources,
        CPPPATH = '..',
        CFLAGS = cflags,
        CXXFLAGS = cxxflags,	
        LIBS = libs)
