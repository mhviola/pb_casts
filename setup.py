"""
Setup configuration for pb_casts package.
"""
from setuptools import setup, find_packages
from pathlib import Path

# Read requirements
requirements = []
req_file = Path(__file__).parent / 'requirements.txt'
if req_file.exists():
    with open(req_file) as f:
        requirements = [line.strip() for line in f if line.strip() and not line.startswith('#')]

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
    packages=find_packages(),
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
    install_requires=requirements,
    include_package_data=True,
)
