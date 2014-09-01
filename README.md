
Normalmap generator
===

Using this tool you can create the following:

![Diffuse Wall](img/diffuse.jpg) ![Normal wall](img/normal.jpg)

Install & libraries
---

You need image magick core to be installed. Then just do:

    $ waf configure build
    Setting top to                           : /home/usr/proj/normalmap 
    Setting out to                           : /home/usr/proj/normalmap/build 
    Checking for 'gcc' (c compiler)          : /usr/bin/gcc 
    Checking for program pkg-config          : /usr/bin/pkg-config 
    Checking for 'MagickCore'                : yes 
    'configure' finished successfully (0.044s)
    'clean' finished successfully (0.004s)
    Waf: Entering directory `/home/rich/proj/normalmap/build'
    [1/3] c: src/normalmap.c -> build/src/normalmap.c.1.o
    [2/3] c: src/main.c -> build/src/main.c.1.o
    [3/3] cprogram: build/src/normalmap.c.1.o build/src/main.c.1.o -> build/normalmap
    Waf: Leaving directory `/home/rich/proj/normalmap/build'
    'build' finished successfully (0.179s)

The binary can then be copied from `build/normalmap`.

Attribution
---

The algorithm to compute the normal map is taken from the gimp normal map plugin. (https://code.google.com/p/gimp-normalmap/).
I found no simple way to create such a normal map via command line.

License
---

Like the gimp normal map plugin GPLv2
