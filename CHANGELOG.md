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
