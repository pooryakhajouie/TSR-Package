from setuptools import setup, find_packages

setup(
    name='tsr_package',
    version='0.1.1',
    packages=find_packages(),
    install_requires=[
        'pandas',
        'numpy',
        'biopython',
        'joblib',
    ],
    author='Poorya Khajouie',
    author_email='poorya.khajouie1@louisiana.edu',
    description='A package for retrieving PDB files and generating key/triplet files for protein analysis.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/pooryakhajouie/TSR-Package',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)
