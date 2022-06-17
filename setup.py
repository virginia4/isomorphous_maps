from setuptools import setup
import pathlib

setup(
   name='stat_ded',
   version='0.1.0',
   author='Alisia Fadini, Virginia Apostolopoulou',
   author_email='apostol.virginia@gmail.com',
   packages=['stat_ded'],
   package_dir={'stat_ded': 'stat_ded'},
   python_requires=">=3.8, <4",
   # scripts=['bin/script1','bin/script2'],
   # install_requires=["cctbx-base"],
   license='LICENSE.txt',
   description='An awesome package that does something',
   long_description=open('README.txt').read()
)
















