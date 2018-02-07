QuPath Bio-Formats Extension
============================

The **QuPath Bio-Formats Extension** makes it possible to use the popular [**Bio-Formats**](http://www.openmicroscopy.org/site/products/bio-formats) library for reading images in [**QuPath**](http://qupath.github.io).

This not only has the advantage of providing support for a wider range of images - including some 16-bit, multichannel, multidimensional data - but it can also help resolve some troubles that come from using native libraries with Java (e.g. when linking up QuPath and [MATLAB](https://github.com/qupath/qupath-matlab-extension)).


## Installation

First, you will need to download two JAR files:

1. The **Bio-Formats Package** from [here](http://www.openmicroscopy.org/site/products/bio-formats/downloads)
2. The **QuPath Bio-Formats Extension** from [here](https://github.com/petebankhead/qupath-bioformats-extension/releases/latest)

Then simply drag both of them onto the main QuPath window while QuPath is running.  If you have not previously installed any QuPath extensions, you will be prompted to select a directory to store them in (or just use the default).

> If you have older versions that can't be overwritten after dragging the .jar files onto QuPath, you can open your QuPath's extensions directory in Explorer/Finder and copy the files in manually there.


## Performance

Starting with v0.0.4, options have been added to QuPath preference panel that can help improve the performance when reading complex images.

It is **highly recommended** to check these out.  More details are available by hovering the cursor over the options, or on [Google Groups](https://groups.google.com/d/msg/qupath-users/78PpZuu2J1s/su6ZjY0mAgAJ).


## Source code

You should find the source code for the **QuPath Bio-Formats Extension** alongside the extension download.

The source code for Bio-Formats is at https://github.com/openmicroscopy/bioformats