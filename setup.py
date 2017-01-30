from setuptools import setup

setup(name='offsetbasedgraph',
      version='0.1.3',
      description='Offset based graph',
      url='http://github.com/uio-cels/offsetbasedgraph',
      author='Ivar Grytten and Knut Rand',
      author_email='',
      license='MIT',
      packages=['offsetbasedgraph'],
      zip_safe=False,
      install_requires=['pymysql', 'numpy', 'future', 'biopython'])