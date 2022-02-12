from setuptools import setup

setup(name='radial-distance-layout',
      version='0.1.1',
      description='Generates a radial layout for trees whose nodes are associated with a distance to the root.',
      url='https://bitbucket.org/bfmaier/radial-distance-layout',
      author='Benjamin F. Maier',
      author_email='bfmaier@physik.hu-berlin.de',
      license='MIT',
      packages=['radial_distance_layout'],
      install_requires=[
          'numpy>=1.17',
          'networkx>=2.0',
      ],
      zip_safe=False)
