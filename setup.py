from setuptools import setup

setup(
    name='metadynminer',
    version='0.2.0',
    description="Python package for efficient analysis of HILLS files generated by Plumed metadynamics simulations. ",
    url='https://github.com/Jan8be/metadynminer.py',
    author='Jan Beránek',
    author_email='Jan1.Beranek@vscht.cz',
    license='GPL-3.0',
    packages=['metadynminer'],
    install_requires=['numpy>=1.21.6',
                      'matplotlib>=3.5.3',
                      'pandas>=1.3.5',
                      'pyvista>=0.38.5'
                      ],

    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',  
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
)
