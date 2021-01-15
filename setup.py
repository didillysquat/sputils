try:
    from setuptools_conda import dist_conda
    cmdclass = {'dist_conda': dist_conda}
except ImportError:
    cmdclass = {}
import setuptools
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="sputils",
    version='0.0.1',
    author="Benjamin C C Hume",
    author_email="didillysquat@gmail.com",
    description="Package including modules to work with SymPortal outputs.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/didillysquat/sputils",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: MacOS",
        "Operating System :: POSIX :: Linux"
    ],
    license='GPL-3.0-only',
    python_requires='>=3.6',
    install_requires=['pandas', 'matplotlib', 'numpy'],
    scripts=['scripts/spbars.py'],
    project_urls={
        'Bug Reports': 'https://github.com/didillysquat/sputils/issues',
        'Source': 'https://github.com/didillysquat/sputils',
    }
)
