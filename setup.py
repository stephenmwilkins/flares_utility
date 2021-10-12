from setuptools import setup, find_packages


# with open('README.rst') as f:
#     readme = f.read()
#
# with open('LICENSE') as f:
#     license = f.read()

setup(
    name='flares_utility',
    version='0.1.0',
    description='for analysing FLARES simulations',
    # long_description=readme,
    author='Stephen Wilkins',
    author_email='s.wilkins@sussex.ac.uk',
    url='https://github.com/stephenmwilkins/flares_utility',
    # license=license,
    packages=find_packages(exclude=('examples', 'useful_plots', 'master_file_scripts')),
    package_data={'flares_utility': ['data/*']}
)
