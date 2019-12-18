# The Biofysica Toolbox

## Installing the Biofysica Toolbox

### Prerequisites

**A working git command line client**

If you don't have one:

On Windows systems get git from [Git-Scm](https://git-scm.com/download/win)

On MacOS systems get git from [Git-Scm](https://git-scm.com/download/mac)
(or use [Homebrew](https://brew.sh), [MacPorts](https://www.macports.org) or the Command Line Developer tools by running xcode-select
```
$ xcode-select --install
```

On Unix/Linux systems use your package manager to install it.

**A working Python 3 environment.**

You can get this from www.python.org

### Installing

Download the [DownloadBiofysicaToolbox.m](https://gitlab.science.ru.nl/marcw/biofysica/blob/master/DownloadBiofysicaToolbox.m) script.
https://gitlab.science.ru.nl/marcw/biofysica/blob/master/DownloadBiofysicaToolbox.m

Run DownloadBiofysicaToolbox from Matlab and follow the instructions, e.g.
```
DownloadBiofysicaToolbox  % install into the default location 
```
or
```
DownloadBiofysicaToolbox.m('/Users/YOU/Documents/MATLAB/myspeciallocation')
```

Add the following two lines to your startup.m file:
```
addpath('/where/your/toolbox/was/installed');
init_biofysica;
```

If there are any lines referring to old installations in startup.m remove them now.

Never manually add the whole toolbox tree to your Matlab path. Neither from the
menu, nor using addpath(genpath(...)) constructs. init_biofysica wil take care of this and
keep everything up to date between updates.


# Updating the Biofysica Toolbox

Run UpdateBiofysicaToolbox from Matlab




