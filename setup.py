# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 17:30:33 2020

Github: https://github.com/tjczec01

@author: Travis J Czechorski

E-mail: tjczec01@gmail.com

"""

import os
from setuptools import setup, find_packages

cwd = os.getcwd()
dir_path = os.path.dirname(os.path.realpath(__file__))

description_chemsys = str("""Calculates the butcher tableau for a given order.""")


with open(r"{}\README.md".format(cwd), "r") as fh:
    long_description = fh.read()

# with open(r'{}\LICENSE.txt'.format(cwd)) as f:
#     license = f.read()


setup(
    name="butchertableau", 
    version="1.0.6",
    author="Travis Czechorski",
    author_email="tjczec01@gmail.com",
    description="{}".format(description_chemsys),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url=r"https://github.com/tjczec01/butcher",
    packages = find_packages(),
    keywords = ['chemical engineering', 'chemistry', 'engineering',
                'chemical reactions', 'Butcher Tableau', "Radau", "Numerical"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent"],
    python_requires='>=3.6')
