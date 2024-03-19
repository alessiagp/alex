from setuptools import setup

setup(
   name='alex',
   version='1.0',
   description='Postprocessing and analysis of data obtained from the EXCOGITO software',
   author='Alessia Guadagnin Pattaro',
   author_email='alessia.guadagnin@unitn.it',
   packages=['alex'],  #same as name
   install_requires=['numpy', 'pandas', 'scipy', 'mdtraj'], #external packages as dependencies
)
