import sys,os
import time
from sys import stdout,stderr
from glob import glob
import platform

from distutils.core import setup,Extension

import distutils.sysconfig
import os

main_libdir=distutils.sysconfig.get_python_lib()
pylib_install_subdir = main_libdir.replace(distutils.sysconfig.PREFIX+os.sep,'')
pylib_install_subdir = pylib_install_subdir.replace('dist-packages','site-packages')

if not os.path.exists('ups'):
    os.mkdir('ups')
tablefile=open('ups/erutil.table','w')
tab="""
# The default version of this file will be overwritten on setup to include
# paths determined from the python version.  This is useful to have in place
# though so that dependencies can be checked *before* installation.  Currently
# there are no required dependencies, so this is somewhat moot.

setupOptional("python")
envPrepend(PYTHONPATH,${PRODUCT_DIR}/%s)
""" % pylib_install_subdir 
tablefile.write(tab)
tablefile.close()



packages = ['erutil','erutil.des','erutil.healpix']
ext_modules = []
include_dirs=[]

#try:
#    import numpy
#    include_dirs=[numpy.get_include()]
#    include_dirs += ['erutil/include']
#    have_numpy=True
#except:
#    have_numpy=False
#    ext_modules=[]
#    include_dirs=[]
#
#    stdout.write('Numpy not found:  Not building C extensions\n')
#    time.sleep(5)
#
#if platform.system()=='Darwin':
#    extra_compile_args=['-arch','i386','-arch','x86_64']
#    extra_link_args=['-arch','i386','-arch','x86_64']
#else:
#    extra_compile_args=[]
#    extra_link_args=[]

try:
    import esutil
except:
    stdout.write('esutil not found: please install from http://code.google.com/p/esutil/')
    sys.exit(2)

exec(open('erutil/version.py').read())

# data_files copies the ups/esutil.table into prefix/ups
setup(name='erutil',
      description="Eli Rykoff's Python Utilities",
      url='https://github.com/erykoff/erutil',
      packages=packages,
      version=__version__,
      data_files=[('ups',['ups/erutil.table'])],
      ext_modules=ext_modules,
      include_dirs=include_dirs)
