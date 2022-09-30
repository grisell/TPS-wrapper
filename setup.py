from distutils.core import setup
setup(
    name = "pygromacstps",
    packages = ["pygromacstps","pygromacstps.orderparameter"],
    package_dir = {"pygromacstps": "pygromacstps"},
    version = "0.0.1",
    description = "Python Gromacs Path Sampling",
    author = "Wolfgang Lechner",
    author_email = "wolfgang.lechner@gmail.com",
    url = "http://homepage.univie.ac.at/wolfgang.lechner",
    download_url = "",
    keywords = ["gromacs", "chemistry", "physics", "scientific"],
    classifiers = [
        "Programming Language :: Python",
        "Development Status :: 4 - Beta",
        "Environment :: Other Environment",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering"
        ],
    long_description = """\
    
    """
)
