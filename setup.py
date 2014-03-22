from setuptools import setup

setup(name='pytides',
      description='Tidal analysis and prediction library.',
      version='0.0.4',
      author='Sam Cox',
      author_email='sam.cox@cantab.net',
      url='http://github.com/sam-cox/pytides',
      packages=['pytides'],
      install_requires=['numpy>=1.8','scipy>=0.11'],
      license='MIT')
