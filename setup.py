import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    dependencies = fh.readlines()

setuptools.setup(
    name="QMCblip",
    version="1.0.5",
    author="Emiel Slootman",
    author_email="e.slootman@esciencecenter.nl",
    description="A small package to couple Quantum Monte Carlo codes to Machine Learning Force Fields through ASE.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/NLESC-JCER/QMCblip",
    project_urls={
        "Bug Tracker": "https://github.com/NLESC-JCER/QMCblip/issues",
    },
    license="Apache Software License 2.0",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
    test_suite='tests',
    packages=setuptools.find_packages(),
    install_requires=dependencies,
    python_requires=">=3.7, <3.10",
)