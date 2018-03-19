from setuptools import setup

setup(name='offsetbasedgraph',
      version='2.0.10',
      description='Offset based graph',
      url='http://github.com/uio-cels/offsetbasedgraph',
      author='Ivar Grytten and Knut Rand',
      author_email='',
      license='MIT',
      packages=['offsetbasedgraph'],
      zip_safe=False,
      install_requires=['pymysql', 'numpy', 'future', 'scipy', 'pyvg'],
      classifiers=[
            'Programming Language :: Python :: 3'
      ]

      )


"""
To update package:
#Update version number manually in this file

sudo python3 setup.py sdist
sudo python3 setup.py bdist_wheel
twine upload dist/offsetbasedgraph-2.0.10.tar.gz
"""