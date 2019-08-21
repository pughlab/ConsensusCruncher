from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='ConsensusCruncher',
      version='0.0.1',
      author="Nina Wang",
      author_email="nina.tt.wang@gmail.com",
      description="Barcoding-based error suppression algorithm",
      long_description=long_description,
      long_description_content_type="text/markdown",
      url="https://github.com/pughlab/ConsensusCruncher",
      packages=["consensuscruncher"],
      install_requires=['numpy', 'pandas', 'pysam', 'Biopython', 'matplotlib'],
      classifiers="Programming Language :: Python :: 3"
      )
