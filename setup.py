import os
from setuptools import setup, find_packages



def read(fname):
    with open(os.path.join(os.path.dirname(__file__), fname)) as f:
        out = f.read()
    return out


def update_version():
    ver = read('VERSION').rstrip('\n')
    fdn = os.path.join(os.path.dirname(__file__), 'semiannotate')
    fn = os.path.join(fdn, 'version.py')
    with open(fn, 'wt') as f:
        f.write('version = "{:}"'.format(ver))
    return ver


setup(
    name="semiannotate",
    version=update_version(),
    author="Fabio Zanini",
    author_email="fabio.zanini@fastmail.fm",
    description="Semi-supervised k-nearest neighbor graphs.",
    license="MIT",
    keywords="graph semi-supervised",
    url="https://github.com/iosonofabio/semiannotate",
    packages=['semiannotate'] + ['semiannotate.' + s for s in find_packages(where='semiannotate')],
    long_description='''
    Atlas-based cell type annotation, with freedom to be queer.


    NOTE: The module leidenalg to perform graph-based clstering is released
    under the GLP3 license. You agree with those licensing terms if you use
    leidenalg within SemiAnnotate.
    ''',
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
        "License :: OSI Approved :: MIT License",
    ],
)
