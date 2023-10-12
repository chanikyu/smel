## setup.py
from setuptools import find_packages, setup

with open("README.md", "r") as fh:
    long_description = fh.read()

NAME = 'smel'
VERSION = '0.0.3'
DESCRIPTION = 'Automatic analysis pipeline'
AUTHOR = 'Kyu-Chan Lee'
EMAIL = 'lkcmsd@gmail.com'
GIT = 'https://github.com/KyuKyu-dplearn/smel'
LICENSE = 'KyuKyu'

INSTALL_REQUIRES = [
    'biopython==1.81',
    'conorm==1.2.0',
    'pyfiglet==0.8.post1',
    'termcolor==2.3.0',
    'wget==3.2.0',
]

PACKAGES = find_packages()

ENTRY_POINTS = {
    'console_scripts': [
        'smel = bin.smel:main',
    ],
}

setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    author=AUTHOR,
    author_email=EMAIL,
    url=GIT,
    install_requires=INSTALL_REQUIRES,
    packages=PACKAGES,
    entry_points=ENTRY_POINTS,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
)
