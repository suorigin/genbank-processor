from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith('#')]

setup(
    name="genbank-processor",
    version="1.2.0",  # 更新版本号
    author="suorigin",
    author_email="1436636379@qq.com",
    description="A comprehensive tool for processing genomic data from NCBI with CDS/protein extraction capabilities",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/genbank-processor",
    packages=find_packages(include=['genbank_processor', 'genbank_processor.*']),
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
    python_requires=">=3.8",
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "genbank-processor=genbank_processor.main:main",
        ],
    },
    include_package_data=True,
)
