## Version 0.0.7

* New calculation of downsample factors
* Support for 32-bit floating point images
* 8-bit indexed color images no longer fail, but treated as regular 8-bit (necessary for some .lif files)
* Revised parsing of z/t values
* New `dumpMetadata()` method to access OME-XML to help with debugging
* New `getChannelName(channel)` method
* Introduced `DummyMetadata()` when creating secondary readers for multithreading
* Removed autoscaling with use of non-RGB `getBufferedThumbnail`
* Improved parsing of objective magnification
* Added new code to test QuPath's use of Bio-Formats
* Added Bio-Formats version number under 'Help &rarr; Installed extensions'
* An error is shown when starting QuPath is 'bioformats_package.jar' is missing

## Version 0.0.6

* Switch to Gradle
* Added workaround for VSI reading bug that meant channels & z-slices were sometimes confused (available in preferences, turned off by default)
* Added fix where 4-channel 8-bit images could be converted to RGB, losing the 4th channel information for some commands
* Added fix where some formats (e.g..ims) would have 8-bit RGB data, but this was treated as multichannel; now such cases are treated as packed (A)RGB, giving better performance & more consistent behavior


## Version 0.0.5

* Fixed bug that prevented handling images with more than 4 channels
* Refined preferences to give more control over performance-related options


## Version 0.0.4

This was a major revision incorporating many bug-fixes and performance improvements, including:
* Better support for a range of image types, especially CZI, NDPIS and VSI
* Bio-Formats' memoization can now be used to improve performance when (re)opening images
* Improved parallelization when reading image tiles - should be faster, require less memory, and turn itself off when it would do more harm than good
* A new set of Bio-Formats options have been added to the QuPath preference panel, including the ability to always use/skip Bio-Formats for specific file formats


## Version 0.0.3

* Fixed bug when resizing image regions


## Version 0.0.2

* Fixed bugs that prevented reading non-RGB whole slide images
* Added nearest-neighbor interpolation method to fix accidental intensity rescaling if resized image regions were requested


## Version 0.0.1

* First available version.
