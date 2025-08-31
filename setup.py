from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = fh.read().splitlines()

setup(
    name="genbank-processor",
    version="1.0.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="A tool for processing, filtering and downloading genomic data from NCBI",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/genbank-processor",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",  
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10", 
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.8",  # 与您的测试矩阵保持一致
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "genbank-processor=genbank_processor.main:main",
        ],
    },
    include_package_data=True,
)
