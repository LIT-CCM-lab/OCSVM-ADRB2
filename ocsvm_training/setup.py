from setuptools import setup

setup(name='ocsvm_training',
      version='0.1',
      description='Package to estimate the nu parameter for OCSVM classifiers using a trained kernel',
      url='https://github.com/LIT-CCM-lab/OCSVM-ADRB2',
      author='Luca Chiesa',
      author_email='luca.chiesa@unistra.com',
      license='MIT',
      packages=['ocsvm_training'],
      install_requires=['numpy', 'matplotlib', 'kneed', 'scipy'],
      include_package_data=True,
      zip_safe=False)