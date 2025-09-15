from setuptools import setup, find_packages

setup(
    name='SpaMeanImpute',
    version='0.1.0',
    description='Spatial Transcriptomics Imputation Tool (SpaMean-Impute)',
    long_description=open('README.md', encoding='utf-8').read(),
    long_description_content_type='text/markdown',
    author='Fahim Hafiz',
    author_email='fahimhafiz@cse.uiu.ac.bd',
    url='https://github.com/FahimHafiz/SpaMean-Impute',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],  
    packages=find_packages(),
    install_requires=[
        'scanpy',
        'scipy',
        'numpy',
        'scikit-learn',
        'tqdm',
        'pandas'
    ],
    python_requires='>=3.7',
    entry_points={
        'console_scripts': [
            'spamean-impute=spa_mean_impute.cli:main',
        ],
    },
)