import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="visual-peaks-AP",
    version="0.0.1",
    author="Aleksandr Fedorov",
    author_email="anfedorov@edu.hse.ru",
    description="Package to compute average precision (AP) for visual peaks in the ChIP-seq experiments.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/alnfedorov/visual-peaks-AP",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering"
    ],
    python_requires='>=3.6',
)
