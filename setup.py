from setuptools import setup

setup(
    name="breakinator",
    version="1.0.0",
    author="Jakob Heinz",
    author_email="jheinz@g.harvard.edu",
    description="Detection of foldback and chimeric read artifacts in PAF files",
    url="https://github.com/jheinz27/breakinator",
    py_modules=["breakinator"], 
    python_requires='>=3.7', 
    entry_points={'console_scripts': ["breakinator=breakinator:main"], },
    )
