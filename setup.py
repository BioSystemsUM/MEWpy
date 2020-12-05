from setuptools import setup, find_packages

files = ["model/data/*"]

with open('README.rst') as readme_file:
    readme = readme_file.read()

requirements = ['optlang<1.4.6', 'python-libsbml', 'inspyred', 'reframed',
                'cobra', 'jmetalpy', 'cobamp', 'networkx']

test_requirements = requirements

setup(
    name='mewpy',
    version='0.1.2',
    python_requires='>=3.6',
    package_dir={'': 'src'},
    packages=find_packages('src'),
    package_data={"": ["*.xml", "*.csv", "*.txt"], 'mewpy': files},
    include_package_data=True,
    zip_safe=False,
    install_requires=requirements,
    author='BiSBII CEB University of Minho',
    author_email='vpereira@ceb.uminho.pt',
    description='mewpy - Metabolic Engineering in Python ',
    license='Apache License Version 2.0',
    keywords='strain optimization',
    url='https://github.com/BioSystemsUM/mewpy/',
    long_description=open('README.rst').read(),
    test_suite='tests',
    tests_require=test_requirements,
    classifiers=[
        'Topic :: Utilities',
        'Programming Language :: Python :: 3.7',
    ],
)
