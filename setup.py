import os, platform, pkg_resources
from distutils.core import setup, Extension
from Cython.Distutils import build_ext

import autowrap;

VERSION = (0, 0, 2);

with open("README.md", "r") as fh:
    long_description = fh.read()

data_dir = pkg_resources.resource_filename("autowrap", "data_files")
include_dir = os.path.join(data_dir, "autowrap")

# if(not os.path.exists(os.path.abspath('./py_chenhancc.cpp'))):
#     import subprocess;
#     source_files_dir = os.path.abspath("./src");
#     source_file_pxd = os.path.join(source_files_dir, "chenhancc.pxd");
#     out_file_pyx = os.path.join(source_files_dir, "py_chenhancc.pyx");
#     
#     subprocess.Popen(" ".join(["autowrap", "--out", out_file_pyx, source_file_pxd]), cwd=os.path.abspath("./src"));
# #     autowrap --out py_chenhancc.pyx chenhancc.pxd

ext = Extension("py_goicp",
                sources = ['src/py_goicp.cpp'],
                language="c++",
                extra_compile_args=["-std=c++14"],
                extra_link_args=["-std=c++14"],
                include_dirs = [include_dir, data_dir],
               )

setup(cmdclass={'build_ext':build_ext},
      name="py_goicp",
      version="%d.%d.%d" % VERSION,
      ext_modules = [ext],
      install_requires=[
          'autowrap',
      ],
      author='#0K Srinivasan Ramachandran',
      author_email='ashok.srinivasan2002@gmail.com',
      url='https://github.com/aalavandhaann/goicp_cython',
      maintainer="#0K Srinivasan Ramachandran",
      maintainer_email="ashok.srinivasan2002@gmail.com",
      platforms=["any"],
      description='GO-ICP compiled using Cython to use in python',
      license='LICENSE.txt',
      keywords='icp go-icp registration alignment rigid-align rigid-alignment',
      python_requires='>=2',
      long_description=long_description,
      long_description_content_type="text/markdown",
      zip_safe=False,
     )

###AUTOWRAP
#autowrap --out py_chenhan.pyx chenhancc.pxd
#python setup.py build_ext --inplace