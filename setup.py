import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r") as fh:
    dependencies = fh.readlines()

setuptools.setup(
    name="CHAMP-ASE-FLARE",
    version="0.0.1",
    author="Emiel Slootman",
    author_email="e.slootman@esciencecenter.nl",
    description="A small package to couple CHAMP to ASE and FLARE",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/NLESC-JCER/CHAMP-ASE-FLARE",
    project_urls={
        "Bug Tracker": "https://github.com/NLESC-JCER/CHAMP-ASE-FLARE/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Physics",
        "Operating System :: OS Independent",
    ],
    packages=setuptools.find_packages(),
    install_requires=dependencies,
    python_requires=">=3.9",
)