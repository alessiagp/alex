from setuptools import setup, find_packages

setup(
   name='alex',
   version='1.0',
   description='Postprocessing and analysis of data obtained from the EXCOGITO software',
   long_description=open("README.md").read() if "README.md" else "",
   author='Alessia Guadagnin Pattaro',
   author_email='alessia.guadagnin@unitn.it',
   packages=find_packages(),  #same as name
   entry_points = {
    'console_scripts': ['alex=alex.run:main']
    },

   install_requires=['numpy', 'pandas', 'scipy', 'mdtraj', 'logging', 'random', 'MDAnalysis'], #external packages as dependencies
   python_requires=">=3.12",
)
