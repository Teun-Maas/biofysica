# The Biofysica Toolbox

## Installing the Biofysica Toolbox

### Prerequisites

**A working git command line client**

If you don't have one:

On Windows systems get git from [Git-Scm](https://git-scm.com/download/win)
Select to INSTALL FOR ALL USERS, ADD IT TO THE PATH, OVERRIDE PATH LENGHT

On MacOS systems get git from [Git-Scm](https://git-scm.com/download/mac)
(or use [Homebrew](https://brew.sh), [MacPorts](https://www.macports.org) or the Command Line Developer tools by running xcode-select
```
$ xcode-select --install
```

On Unix/Linux systems use your package manager to install it.

**A working Python 3 environment.**

You can get this from www.python.org and make sure you install a version that is
compatible with your Matlab version. Google for "matlab python compatibility" to find out which one. Currently (2021-08) Matlab 2021a supports Python 3.7 and 3.8.

Install these python3 libraries:
```
pip3 install  msgpack numpy scipy
```

MacOS: install homebrew from brew.sh
The installation of the libraries is slightly different:
```
brew install python3 numpy scipy
pip3 install msgpack
```

Find out where the currently installed python3 is located.
Unix/MacOS
```
bash:~$ type -a python3
python3 is /usr/local/bin/python3
python3 is /usr/local/bin/python3
python3 is /usr/bin/python3

```

Windows: open a new CMD.EXE window and type
```
C:\> where python
C:\Program Files\Python38\python.exe
C:\Users\you\Appdata.......\python.exe

```
Usually the upmost line is the location of the python program you want.

If you don't have a startup.m in your Matlab [userpath folder](https://www.mathworks.com/help/matlab/ref/userpath.html), create one.

Add the location of python you found above to the TOP of your startup.m, e.g.

```
pyversion('C:\Program Files\Python38\python.exe');
```

### Installing

Open a terminal window with bash in it. (On Windows look for Git bash in the Start menu)
Navigate to your Matlab userdir, e.g. in Windows and clone the git repository containing the biofysica toolbox:

```
$ cd ~/Documents/MATLAB
$ git clone https://gitlab.science.ru.nl/marcw/biofysica.git
```
The toolbox directories are now in ~/Documents/MATLAB/biofysica

Add the following lines to your startup.m file
Windows:
```
addpath('C:\Users\MyName\Documents\MATLAB\biofysica');
init_biofysica;
```
MacOS/Linux:
```
addpath('~/MATLAB/biofysica');
init_biofysica;
```

If there are any lines referring to old installations in startup.m remove them now.

***IMPORTANT:***
Never manually add the whole toolbox tree to your Matlab path. Neither from the
menu, nor using addpath(genpath(...)) constructs. init_biofysica wil take care of this and
keep everything up to date between updates.


Start a fresh Matlab to get going!

### Updating

Usually it is enough to open a (git) bash command window, navigate into the biofysica directory and type the command
```
git pull
```
This will update the contents of the directory tree. After restarting MATLAB everythin shoud be okay.



