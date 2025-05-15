import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
     name='litraj',  
     version='0.1',
     py_modules = ["litraj"],
     install_requires = [
                            "requests",
                            "pandas",
                            "numpy",
                            "scipy",
                            "ase",
                            "tqdm",
                         ],
     author="Artem Dembitskiy",
     author_email="art.dembitskiy@gmail.com",
     description="LiTraj: Li-ion migration dataset",
     key_words = ["Li-ion", "migration-barrier", "DFT", "machine-learning", "dataset", "benchmark", "universal-potential"],
     long_description=long_description,
     long_description_content_type="text/markdown",
     url="https://github.com/AIRI-Institute/LiTraj",
     package_data={"litraj": ["*.txt", "*.rst", '*.md'], 
                    },
     classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: MIT License",
         "Operating System :: OS Independent",
     ],
    include_package_data=True,
    packages=setuptools.find_packages(),
 )