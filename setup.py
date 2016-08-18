from setuptools import setup
import glob
import os.path as op
from os import listdir

__version__ = '0.0.1'

setup_args = {
    'name': 'katalogss',
    'author': 'Patti',
    'license': 'BSD',
    'package_dir': {'katalogss': 'katalogss'},
    'packages': ['katalogss'],
    'scripts': glob.glob('scripts/*'),
    'version': __version__,
    'package_data': {'katalogss': [f for f in listdir('./katalogss/data') if op.isfile(op.join('./katalogss/data', f))]},
    'install_requires': ['numpy', 'scipy', 'astropy>=1.2','scikit-learn','pandas','pidly']
}

if __name__ == '__main__':
    apply(setup, (), setup_args)
