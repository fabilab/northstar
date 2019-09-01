import os
from setuptools import setup, find_packages



def read(fname):
    with open(os.path.join(os.path.dirname(__file__), fname)) as f:
        out = f.read()
    return out


def update_version():
    ver = read('northstar/_version.py').rstrip('\n').split()[-1].strip('"')
    return ver


setup(
    name="northstar",
    version=update_version(),
    author="Fabio Zanini",
    author_email="fabio.zanini@fastmail.fm",
    description="Cell type annotation guided by cell atlases, with freedom to be queer.",
    license="MIT",
    keywords="graph semi-supervised",
    url="https://github.com/iosonofabio/northstar",
    packages=['northstar'] + ['northstar.' + s for s in find_packages(where='northstar')],
    long_description='''
    Atlas-based cell type annotation, with freedom to be queer.


    NOTE: The module leidenalg to perform graph-based clstering is released
    under the GLP3 license. You agree with those licensing terms if you use
    leidenalg within northstar.
    ''',
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
        "License :: OSI Approved :: MIT License",
    ],
)
