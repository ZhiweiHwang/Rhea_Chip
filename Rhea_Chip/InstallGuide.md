#Objective
Install all software packages required to follow the GATK Best Practices.

#Prerequisites
To follow these instructions, you will need to have a basic understanding of the meaning of the following words and command-line operations. If you are unfamiliar with any of the following, you should consult a more experienced colleague or your systems administrator if you have one. There are also many good online tutorials you can use to learn the necessary notions.

> * Basic Unix environment commands
> * Binary / Executable
> * Compiling a binary
> * Adding a binary to your path
> * Command-line shell, terminal or console

#Software library
You will also need to have access to an ANSI compliant C++ compiler and the tools needed for normal compilations (make, shell, the standard library, tar, gunzip). These tools are usually pre-installed on Linux/Unix systems. **On MacOS X, you may need to install the MacOS Xcode tools. See [https://developer.apple.com/xcode/](https://developer.apple.com/xcode/) for relevant information and software downloads. **The XCode tools are free but an AppleID may be required to download them.

Rhea_Chip requires python (version >= 2.7) and Java (version >= 1.7). All Linux/Unix and MacOS X systems should have a JRE pre-installed, but the version may vary. To test your Java version, run the following command in the shell:

```
python -V
java -version
```
 
This should return a message along the lines of ”java version 1.7.0_25” as well as some details on the Runtime Environment (JRE) and Virtual Machine (VM). If you have a version other than 1.7.x, be aware that you may run into trouble with some of the more advanced features of the Picard and GATK tools. The simplest solution is to install an additional JRE and specify which you want to use at the command-line. To find out how to do so, you should seek help from your systems administrator.

#Software packages
> * BWA
> * SAMtools
> * Picard
> * Genome Analysis Toolkit (GATK)
> * SOAPnumke
> * snpEff
> * python libraries

*Note that the version numbers of packages you download may be different than shown in the instructions below. If so, please adapt the number accordingly in the commands.*

Installation

```
export $Rhea_Chip_Home=$HOME
python setup install
python Rhea_Chip/Rhea_Chip/install.py
```


## Important note

We do not provide all database files which will be used in our analysis, everyone can build the use databases refer on our demo.
RheaChip is Open-source, but some softwares in our pipeline may be commercial use, THIS IS WHAT WE NEED TO PAY ATTENTION !!

