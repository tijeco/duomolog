from setuptools import setup

setup(name='duomolog',
      version='0.1',
      description='simples example in the python world',
      url='https://github.com/gdamjan/hello-world-python-package',
      author='gdamjan',
      author_email='gdamjan',
      license='MIT',
      packages=['duomolog'],
      zip_safe=False,
      entry_points = {
          'console_scripts': ['duomolog=duomolog.duomolog:main'],
      }
)