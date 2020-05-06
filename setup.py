from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="wq_toolbox",
    version="0.0.1",
    description="PySWMM toolbox for modelling any pollutant generation or removal method",
    author="Brooke Mason, Abhiram Mullapudi",
    author_email="bemason@umich.edu, abhiramm@umich.edu",
    long_description=long_description,
    url="https://github.com/kLabUM/wq_toolbox",
    packages=['wq_toolbox'],
    package_data=setuptools.find_packages(),
    classifiers=[
        "Programming Lanageuage :: Python :: 3",
        "License :: OSI Apprpoved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "numpy",
        "pyswmm",
        "scipy",
    ],
    python_requires='>=3.6',
    )

