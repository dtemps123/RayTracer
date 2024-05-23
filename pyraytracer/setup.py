from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='pyraytracer',
      version='0.1',
      description='Optical ray tracing code that propagates rays through a surface based geometry.',
	    long_description=long_description,
      long_description_content_type="text/markdown",
      url='https://github.com/dtemps123/RayTracer',
      author='C. Eric Dahl, Dylan Temples',
      author_email='dtemples@fnal.gov',
      license='MIT',
      packages=setuptools.find_packages(),
      classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
      ],
      python_requires='>=3.6',
      install_requires=[
        'h5py',
        'numpy',
        'scipy',
        'pandas',
        'matplotlib',
      ],
      zip_safe=False)
