from setuptools import setup, find_packages

exec(open('fieldpathogenomics/version.py').read())

setup(
    name = "FieldPathogenomics",
    version = __version__,
    author = "Daniel Bunintg",
    author_email = "daniel.bunting@earlham.ac.uk",
    url = "https://github.com/dnlbunting/FieldPathogenomics/",
    packages=find_packages(),

)
