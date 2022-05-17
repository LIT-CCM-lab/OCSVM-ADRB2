from setuptools import setup

setup(name='pyichem',
      version='0.1',
      description='Package to interact with IChem from a python interface',
      url='https://github.com/LIT-CCM-lab/OCSVM-ADRB2',
      author='Luca Chiesa',
      author_email='luca.chiesa@unistra.com',
      license='MIT',
      packages=['pyichem'],
      install_requires=['numpy', 'pyaml', 'pandas', 'grakel', 'scipy', 'biopandas'],
      include_package_data=True,
      zip_safe=False)