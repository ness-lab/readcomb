import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="readcomb-aays",
    version="0.0.1",
    author="Jiyu Liu and Ahmed Hasan",
    # author_email="author@example.com",
    description="Fast detection of recombinant reads in BAMs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ness-lab/readcomb",
    license="GNU GPLv3+",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Development Status :: 3 - Alpha",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.5',
)
