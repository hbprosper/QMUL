# QMUL
Assorted statistics and ML examples to go with the talk at https://indico.ph.qmul.ac.uk/indico/conferenceDisplay.py?confId=747

## Introduction
The jupyter notebooks here are tutorials that illustrate some of the ideas discussed in the talk given both at PHYSTAT2020 as well as ta Queen Mary University, London (Oct. 20, 2020).

## Installation
The machine learning example can be run on Google's Colab in order to take advantage of their GPUs. But, it can also be run on your local machine. But to do so, you will need to install a set of Python packages. The simplest way to do so is first to install miniconda (a slim version of Anaconda) on your laptop by following the instructions at:

https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html

I recommend installing miniconda3, which comes pre-packaged with Python 3.

Software release systems such as Anaconda (__conda__ for short) make it possible to have several separate self-consistent named *environments* on a single machine, say your laptop. For example, you may need to use Python 2.7.14 sometimes and Python 3.6.7 at other times. If you install software without using *environments* there is the very real danger that the software on your laptop will become inconsistent. Anaconda (and its lightweight companion miniconda) provide a way, for example, to create a software *environment* consistent with Python 2.7.14 and another that is consistent with Python 3.6.7. 

After installing miniconda3, It is a good idea to update conda using the command
```bash
conda update conda
```
Assuming conda is properly installed and initialized on your laptop, you can create an environment, here we call it *python3*, containing the __root__ package from CERN, plus a large subset of the packages in the conda system, using the command>
```bash
conda create -c conda-forge --name python3 root
```
Before pressing __y__ to continue with the installation, scan through the list of packages and identify which of the above are in the list. That way, you will know which ones are missing and need to be installed using the conda install command. For example, as of this writing __pytorch__ is not available by default. In order to install a missing packages, first be sure to choose the conda environment into which the package is to be installed. First activate the desired environment, by doing, for example,
```bash
conda activate python3
```
Later, in order to update root together with a consistent set of packages do
```bash
conda update root
```
taking care to do so in the desired conda environment, here __python3__.

### Other Packages

To install pytorch do
```bash
conda install pytorch torchvision -c pytorch
```
You may also wish to install the rather impressive 3D animation package __vpython__,
```bash
conda install vpython -c vpython
```

If all goes well, you will have installed a rather complete set of amazing high quality *absolutely free* software packages on your system that are consistent with Python 3.6.7.

For some quick help on conda see 

https://uoa-eresearch.github.io/eresearch-cookbook/recipe/2014/11/20/conda/


If you still prefer to do everything by hand, follow the instructions at

https://www.scipy.org/install.html

and 

https://jupyter.org/install

## Tutorials


