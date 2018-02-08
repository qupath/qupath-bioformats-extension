QuPath Bio-Formats Extension
============================

The **QuPath Bio-Formats Extension** makes it possible to use the popular [**Bio-Formats**](http://www.openmicroscopy.org/site/products/bio-formats) library for reading images in [**QuPath**](http://qupath.github.io).

This not only has the advantage of providing support for a wider range of images - including some 16-bit, multichannel, multidimensional data - but it can also help resolve some troubles that come from using native libraries with Java (e.g. when linking up QuPath and [MATLAB](https://github.com/qupath/qupath-matlab-extension)).


## Installation

First, you will need to download two JAR files:

1. The **Bio-Formats Package** from [here](http://www.openmicroscopy.org/site/products/bio-formats/downloads)
2. The **QuPath Bio-Formats Extension** from [here](https://github.com/qupath/qupath-bioformats-extension/releases/latest)

Then simply drag both of them onto the main QuPath window while QuPath is running.  If you have not previously installed any QuPath extensions, you will be prompted to select a directory to store them in (or just use the default).

> If you have older versions that can't be overwritten after dragging the .jar files onto QuPath, you can open your QuPath's extensions directory in Explorer/Finder and copy the files in manually there.


## Performance

Starting with v0.0.4, options have been added to QuPath preference panel that can help improve the performance when reading complex images.

It is **highly recommended** to check these out.  More details are available by hovering the cursor over the options, or on [Google Groups](https://groups.google.com/d/msg/qupath-users/78PpZuu2J1s/su6ZjY0mAgAJ).

![Preferences](prefs.jpg)

A brief description:
* **Enable Bio-Formats**
  * If unchecked, Bio-Formats is disabled and QuPath will try only other image readers
* **Enable Bio-Formats tile parallelization**
  * If checked, Bio-Formats will try to read image tiles in parallel.  This can make browsing whole slide images faster, but can also require more memory and processing power.  If QuPath is able to calculate that the memory requirements would be too high for a specific image, then it will automatically turn this setting off for that image.
* **Enable Bio-Formats channel parallelization (experimental)**
  * If checked, Bio-Formats will try to read multiple channels from the same image tile in parallel.  Where tile parallelization isn't possible because of memory trouble, this can help boost performance for multichannel images.  However, *it does not always work!* I have seen it fail for `.lif` files, but succeed (and help) for some `.czi` files.  If it fails QuPath will try to fall back on requesting channels sequentially - so there are no known adverse effects at the time of writing, but use with caution.
* **Bio-Formats memoization time (ms)**
  * Bio-Formats can sometimes *vastly* speed up opening files & read tiles by writing `.bfmemo` cache files alongside the images.  This setting controls how slow opening an image should be before Bio-Formats will write the cache file (in ms).  Set this to a negative value, or an astronomically high value, if you want to turn off writing `.bfmemo` files altogether (e.g. if images are stored on another server, where new files should not be written).
* **Bio-Formats memoization directory**
  * Specify a directory where Bio-Formats should write the `.bfmemo` files, rather than writing them alongside the corresponding images.
* **Always use Bio-Formats for specified image extensions**
  * QuPath support multiple libraries for opening image formats.  In the case where multiple libraries support the same image formats, the choice of which library is used (e.g. OpenSlide, Bio-Formats) can sometimes seem like a lottery.  Enter image extensions here if you always want Bio-Formats to be preferred for that format.
* **Never use Bio-Formats for specified image extensions**  
  * Like the previous option, but always avoid using Bio-Formats for a specific extension.

> **A example of a time when these final two options can matter:**  OpenSlide can generally read `.ndpi` images much faster than Bio-Formats.  However, some `.ndpi` files can contain z-stacks.  If OpenSlide is asked to read a z-stack, then it will only actually return a single plane from the stack - whereas Bio-Formats can read the entire thing.  Therefore if you use `.ndpi` images you *might* want to force QuPath to use one or the other, depending upon whether you anticipate having z-stacks or not.


## Source code

You should find the source code for the **QuPath Bio-Formats Extension** alongside the extension download.

The source code for Bio-Formats is at https://github.com/openmicroscopy/bioformats
