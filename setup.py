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
    url="https://github.com/northstaratlas/northstar",
    packages=['northstar'] + ['northstar.' + s for s in find_packages(where='northstar')],
    long_description='''
    Cell type annotation guided by cell atlases, with freedom to be queer.

    See https://github.com/northstaratlas/northstar for the project's website.
    ''',
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Utilities",
        "License :: OSI Approved :: MIT License",
    ],
)
