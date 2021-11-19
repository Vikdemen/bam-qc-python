from setuptools import setup

setup(
    name='bamqc',
    version='0.1.0',
    package_dir={'': 'src'},
    packages=['bamqc'],
    url='https://github.com/Vikdemen',
    license='MIT',
    author='Demenev Viktor',
    author_email='viktor.demen@gmail.com',
    description='A small tool showcasing parallel reading of BAM files. Calculates average GC content of uniquely '
                'mapped reads.',
    python_requires='>= 3.9',
    entry_points={
        'console_scripts': ['bamqc=bamqc.bamqc:main'],
    }
)
