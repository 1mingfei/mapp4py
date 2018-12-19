from distutils.core import setup, Extension
import os
from numpy import __path__ as npy_path
npy_lib=['npymath']
npy_lib_dir = npy_path[0] + '/core/lib'
npy_include_dir = npy_path[0] + '/core/include'


mpi_lib=['mpi','mpi_cxx','open-rte','open-pal']
mpi_lib_dir = '/sw/arcts/centos7/openmpi/1.10.2-intel-17.0.1-1/lib'
mpi_include_dir = '/sw/arcts/centos7/openmpi/1.10.2-intel-17.0.1-1/include'


os.environ["CC"] = 'icpc'
os.environ["CXX"] = 'icpc'



cpp_files=[]
cpp_files += ['src/'+ each for each in os.listdir('src') if each.endswith('.cpp')]

module1 = Extension('mapp',
                    include_dirs = [npy_include_dir,mpi_include_dir],
                    libraries =npy_lib+mpi_lib,
                    library_dirs = [npy_lib_dir,mpi_lib_dir],
                    sources = cpp_files,
		    extra_compile_args=['-std=c++11','-O3'])

setup (name ='mapp4py',
       version = '0.0.0',
       description = 'MIT Atomistic Parallel Package',
       author = 'Sina Moeini',
       author_email = 'sinam@mit.edu',
       url = 'https://github.com/sinamoeini/mapp4py',
       ext_modules = [module1])

