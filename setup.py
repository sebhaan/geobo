## To install locally: python setup.py build && python setup.py install
## (If there are problems with installation of the documentation, it may be that
##  the egg file is out of sync and will need to be manually deleted - see error message
##  for details of the corrupted zip file. )
##
## To push a version through to pip.
##  - Make sure it installs correctly locally as above
##  - Update the version information in this file
##  - python setup.py sdist upload -r pypitest  # for the test version
##  - python setup.py sdist upload -r pypi      # for the real version
##
## (see https://gist.github.com/davydany/b08acef08f75fe297e13ae4d24ce9f4d)


from setuptools import setup, find_packages
from numpy.distutils.core import setup, Extension
from os import path
import os
import subprocess
import io

## in development set version to none and ...
PYPI_VERSION = '0.1.0'

# Return the git revision as a string (from numpy)
def git_version():
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        out = subprocess.Popen(cmd, stdout = subprocess.PIPE, env=env).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', '--short', 'HEAD'])
        GIT_REVISION = out.strip().decode('ascii')
    except OSError:
        GIT_REVISION = "Unknown"

    return GIT_REVISION


if PYPI_VERSION is None:
    PYPI_VERSION = git_version()


this_directory = path.abspath(path.dirname(__file__))
with io.open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

packages = find_packages()

if __name__ == "__main__":
    setup(name = 'geobo',
          author            = "Sebastian Haan",
          author_email      = "sebastian.haan@sydney.edu.au",
          url               = "https://github.com/sebhaan/geobo/tree/release",
          version           = PYPI_VERSION,
          description       = "Python package for Multi-Objective Bayesian Optimisation and Joint Inversion in Geosciences",
          long_description  = long_description,
          long_description_content_type='text/markdown',
          install_requires  = ['numpy>=1.15.0', 
                              'scipy>=1.1.0', 
                              'matplotlib>=2.2.2',
                              'scikit_image>=0.14.0',
                              'scipy>=1.1.0',
                              'rasterio>=1.0.20',
                              'pandas>=0.23.4',
                              'pyvista>=0.23.0'],
                              #'skimage>=0.0',
                              #'PyYAML>=5.3'],
          python_requires   = '>=3.6',
          setup_requires    = ["pytest-runner", 'webdav'],
          tests_require     = ["pytest", 'webdav'],
          packages          = ['geobo', 'examples'],
          package_data      = {'geobo': ['*.yaml'],
                               'examples': ['testdata/sample/*',
                                        '*.yaml']},
                                          #'example/testsdata/*',
                                          #'example/testsdata/sample/*',
                                          #'example/testsdata/synthetic/*'] },
          include_package_data = False,
          classifiers       = ['Programming Language :: Python :: 3',
                               'Programming Language :: Python :: 3.4',
                               'Programming Language :: Python :: 3.5',
                               'Programming Language :: Python :: 3.6',
                               'Programming Language :: Python :: 3.7'
                               ]
          )