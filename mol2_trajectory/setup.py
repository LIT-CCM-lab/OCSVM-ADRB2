from setuptools import setup

setup(name='mol2_trajectory',
      version='0.1',
      description='Package to convert and manage MD trajectories as collections of mol2 files',
      url='https://github.com/LIT-CCM-lab/OCSVM-ADRB2',
      author='Luca Chiesa',
      author_email='luca.chiesa@unistra.com',
      license='MIT',
      packages=['mol2_trajectory'],
      install_requires=['numpy', 'pyaml', 'pandas', 'pytraj'],
      include_package_data=True,
      package_data = {'conversion_files': ['mol2_trajectory/conversion_file']},
      zip_safe=False)