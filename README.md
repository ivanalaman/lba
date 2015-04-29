### lba

This is the development code of the R package **lba**.
You should use it if you want to contribute to its development:
testing unreleased versions, fixing bugs, writing code, etc.

To download, check and build it do the following in a terminal emulator:

> git clone  git://github.com/ivanalaman/lba.git

> or

> git clone https://ivanalaman@github.com/lba.git

After to clone it, to check, build and install do the following:
> R CMD check lba

> R CMD build lba

> R CMD INSTALL lba_X.X-X.tar.gz

Or, you can install using devtools package as:

> library(devtools)

> install_github('ivanalaman/lba')

The stable version of this package is available at:
http://cran.r-project.org/web/packages/lba/index.html
