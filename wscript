

def options(c):
    c.load("compiler_c")

def configure(c):
    c.load("compiler_c")
    c.check_cfg(package='MagickCore', args='--cflags --libs', uselib_store='MAGICK_CORE')
    c.env.CFLAGS += ['-std=c11']

def build(c):
    c.program(source='src/normalmap.c src/main.c',
        target='normalmap',
        use='MAGICK_CORE',
        lib='m',
        includes='src')
