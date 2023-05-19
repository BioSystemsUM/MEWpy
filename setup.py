from setuptools import setup, find_packages
import sys

files = ["model/data/*"]

with open('README.rst') as readme_file:
    readme = readme_file.read()

requirements = ['cobra', 'inspyred', 'jmetalpy<=1.5.5',
                'reframed', 'networkx', 'matplotlib<=3.5.0',
                'joblib', 'tdqm', 'httpx<=0.23.0']

setup_requirements = requirements + ['pytest-runner']
test_requirements = requirements + ['pytest', 'cplex']
install_requirements = requirements

setup(
    name='mewpy',
    version='0.1.29',
    python_requires='>=3.6',
    package_dir={'': 'src'},
    packages=find_packages('src'),
    package_data={"": ["*.xml", "*.csv", "*.txt"], 'mewpy': files},
    include_package_data=True,
    zip_safe=False,
    install_requires=install_requirements,
    setup_requires=setup_requirements,
    tests_require=test_requirements,
    author='BiSBII CEB University of Minho',
    author_email='vpereira@ceb.uminho.pt',
    description='mewpy - Metabolic Engineering in Python ',
    license='Apache License Version 2.0',
    keywords='strain optimization',
    url='https://github.com/BioSystemsUM/mewpy/',
    long_description=readme,
    test_suite='tests',
)
