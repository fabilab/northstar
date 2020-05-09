import os
from setuptools import setup, find_packages



def read(fname):
    this_directory = os.path.abspath(os.path.dirname(__file__))
    with open(os.path.join(this_directory, fname)) as f:
        out = f.read()
    return out


def update_version():
    ver = read('northstar/_version.py').rstrip('\n').split()[-1].strip('"')
    return ver


def get_readme():
    return read('README.md')


setup(
    name="northstar",
    version=update_version(),
    author="Fabio Zanini",
    author_email="fabio.zanini@fastmail.fm",
    description="Single cell type annotation guided by cell atlases, with freedom to be queer.",
    license="MIT",
    keywords="graph semi-supervised",
    url="https://github.com/northstaratlas/northstar",
    packages=['northstar'] + ['northstar.' + s for s in find_packages(where='northstar')],
    long_description=get_readme(),
    long_description_content_type='text/markdown',
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Utilities",
        "License :: OSI Approved :: MIT License",
    ],
    install_requires=[
        'numpy',
        'pandas',
        'scikit-learn',
        'python-igraph>=0.8.0',
        'leidenalg>=0.8.0',
        'anndata',
    ],
    setup_requires=[
        'numpy',
        'pandas',
        'scikit-learn',
        'python-igraph>=0.8.0',
        'leidenalg>=0.8.0',
        'anndata',
    ],
    extras_require={
        'atlas-fetcher': [
            'requests',
            'loompy',
        ]
    },
)
