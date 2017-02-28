from setuptools import setup

setup(name='offsetbasedgraph',
      version='1.0.0',
      description='Offset based graph',
      url='http://github.com/uio-cels/offsetbasedgraph',
      author='Ivar Grytten and Knut Rand',
      author_email='',
      license='MIT',
      packages=['offsetbasedgraph'],
      zip_safe=False,
      install_requires=['pymysql', 'numpy', 'future', 'biopython'])


"""
To update package:

sudo python3 setup.py sdist
sudo python3 setup.py bdist_wheel
twine upload dist/*
"""