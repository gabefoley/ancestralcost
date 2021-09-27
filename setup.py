from setuptools import setup
import os
import re


def read_version():
    path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'ancestralcost/__init__.py')
    with open(path, 'r') as fh:
        return re.search(r'__version__\s?=\s?[\'"](.+)[\'"]', fh.read()).group(1)


def readme():
    with open('README.md') as f:
        return f.read()


setup(name='ancestralcost',
      version=read_version(),
      description='',
      long_description=readme(),
      long_description_content_type='text/markdown',
      author='Gabe Foley',
      author_email='g.foley@uq.edu.au',
      url='https://github.com/gabefoley/ancestralcost',
      license='GPL3',
      project_urls={
          "Bug Tracker": "https://github.com/gabefoley/ancestralcost/issues",
          "Documentation": "https://github.com/gabefoley/ancestralcost",
          "Source Code": "https://github.com/gabefoley/ancestralcost",
      },
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Natural Language :: English',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: 3.8',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
      ],
      keywords=['ancestral sequence reconstruction', 'phylogenetics', 'ancestors'],
      packages=['ancestralcost'],
      entry_points={
          'console_scripts': [
              'PROJECT_NAME = ancestralcost.__main__:main'
          ]
      },
      install_requires=['pandas', 'numpy', 'ete3', 'pyqt5'],
      python_requires='>=3.6',
      data_files=[("", ["LICENSE"])]
      )
