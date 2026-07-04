"""
Setup configuration for pb_casts package.
"""
from setuptools import setup
from pathlib import Path

# requirements.txt is read here for reference but install_requires below uses
# a pinned list instead so that conda-only deps (e.g. python-ctd) don't break
# a plain `pip install`.
req_file = Path(__file__).parent / 'requirements.txt'

# Read README
readme = ''
readme_file = Path(__file__).parent / 'README.md'
if readme_file.exists():
    with open(readme_file, encoding='utf-8') as f:
        readme = f.read()

setup(
    name='pb_casts',
    version='0.1.0',
    author='Marisa Viola',
    author_email='your.email@wwu.edu',
    description='CTD data processing package for Padilla Bay research',
    long_description=readme,
    long_description_content_type='text/markdown',
    url='https://github.com/yourusername/pb_casts',
    packages=['pb_casts'],
    package_dir={'pb_casts': '.'},
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Oceanography',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
    ],
    python_requires='>=3.8',
    install_requires=[
        'numpy>=1.20.0',
        'pandas>=1.3.0',
        'matplotlib>=3.4.0',
        'gsw>=3.4.0',
        'xarray>=0.19.0',
        'pyarrow>=5.0.0',
        'openpyxl>=3.0.0',
    ],
    include_package_data=True,
)
