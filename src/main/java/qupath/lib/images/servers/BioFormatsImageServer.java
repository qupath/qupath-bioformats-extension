/*-
 * #%L
 * This file is part of a QuPath extension.
 * %%
 * Copyright (C) 2014 - 2016 The Queen's University of Belfast, Northern Ireland
 * Contact: IP Management (ipmanagement@qub.ac.uk)
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */

package qupath.lib.images.servers;

import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.awt.image.WritableRaster;
import java.io.File;
import java.io.IOException;
import java.nio.channels.ClosedByInterruptException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.WeakHashMap;
import java.util.concurrent.TimeUnit;

import javax.imageio.ImageIO;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ome.units.UNITS;
import ome.units.quantity.Length;
import ome.units.unit.Unit;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.ClassList;
import loci.formats.FormatException;
import loci.formats.IFormatReader;
import loci.formats.ImageReader;
import loci.formats.Memoizer;
import loci.formats.gui.AWTImageTools;
import loci.formats.gui.BufferedImageReader;
import loci.formats.in.NDPIReader;
import loci.formats.meta.IMetadata;
import loci.formats.meta.MetadataStore;
import loci.formats.services.OMEXMLService;
import loci.formats.tiff.IFD;
import qupath.lib.awt.common.AwtTools;
import qupath.lib.awt.images.PathBufferedImage;
import qupath.lib.common.ColorTools;
import qupath.lib.common.GeneralTools;
import qupath.lib.images.PathImage;
import qupath.lib.regions.RegionRequest;

/**
 * QuPath ImageServer that uses the Bio-Formats library to read image data.
 * 
 * See http://www.openmicroscopy.org/site/products/bio-formats
 * 
 * @author Pete Bankhead
 *
 */
public class BioFormatsImageServer extends AbstractImageServer<BufferedImage> {
	
	private static final Logger logger = LoggerFactory.getLogger(BioFormatsImageServer.class);
	
	/**
	 * In the event that we're using multithreading, store multiple readers.
	 */
	private Map<Thread, IFormatReader> map = new WeakHashMap<>();
		
	/**
	 * Image names (in lower case) normally associated with 'extra' images, but probably not representing the main image in the file.
	 */
	private static List<String> extraImageNames = Arrays.asList("overview", "label", "thumbnail", "macro");
	
	/**
	 * Original metadata, populated when reading the file.
	 */
	private ImageServerMetadata originalMetadata;
	
	/**
	 * Optional modified metadata, which may be used to override the originalMetadata.
	 */
	private ImageServerMetadata userMetadata;

	/**
	 * Full path used to identify this image; may contain an extra identifier for a contained sub-image
	 */
	private String path;

	/**
	 * Path to the base image file - will be the same as path, unless the path encodes the name of a specific series, in which case this refers to the file without the series included
	 */
	private String filePath;
	
	/**
	 * The main image reader from which to request pixels.
	 */
	private BufferedImageReader reader;
	
	/**
	 * Preferred downsample values, representing each level of an image pyramid.
	 */
	private double[] downsamples;
	
	/**
	 * The image is 8-bit RGB.
	 */
	private boolean isRGB = true;
	
	/**
	 * Bits per pixel.
	 */
	private int bpp = 0;
	
	/**
	 * Representation of the different time points for a time series.
	 * TODO: Test the use of different time points, if this ever becomes a common use of QuPath.
	 */
	private double[] timePoints = null;
	
	/**
	 * A map linking an identifier (image name) to series number for 'full' images.
	 */
	private Map<String, Integer> imageMap = null;
	
	/**
	 * A map linking an identifier (image name) to series number for additional images, e.g. thumbnails or macro images.
	 */
	private Map<String, Integer> associatedImageMap = null;

	/**
	 * List representing the preferred colors to use for each channel, as packed int values.
	 */
	private List<Integer> channelColors = new ArrayList<>();
	
	/**
	 * Delimiter between the file path and any sub-image names
	 */
	private String delimiter = "::";

	/**
	 * Numeric identifier for the image (there might be more than one in the file)
	 */
	private int series = 0;
	
	/**
	 * QuPath-specific options for how the image server should behave, such as using parallelization or memoization.
	 */
	private BioFormatsServerOptions options;
	
	/**
	 * Size of any memoization file, in bytes.
	 * This is useful in any heuristic aiming to turn of parallel reading for especially heavyweight image readers.
	 */
	private long memoizationFileSize = -1L;
	
	/**
	 * Create an ImageServer using the Bio-Formats library for a specified image path.
	 * 
	 * @param path
	 * @throws FormatException
	 * @throws IOException
	 * @throws DependencyException
	 * @throws ServiceException
	 */
	public BioFormatsImageServer(final String path) throws FormatException, IOException, DependencyException, ServiceException {
		this(path, BioFormatsServerOptions.getInstance());
	}
	
	private static Map<String, BufferedImageReader> cachedReaders = new HashMap<>();

	BioFormatsImageServer(final String path, final BioFormatsServerOptions options) throws FormatException, IOException, DependencyException, ServiceException {
		super();
		
		long startTime = System.currentTimeMillis();

		this.options = options;

		// Zip files tend to be slow for Bioformats to parse... so best not try
		if (path.toLowerCase().endsWith(".zip"))
			throw new FormatException("ZIP not supported");

		// Warning, because my interpretation of what Bio-Formats is doing may be wrong...
		// leading to some unfortunate bugs
		logger.warn("QuPath Bio-Formats extension is in beta!  Be watchful for bugs...");
				
		// Create variables for metadata
		int width = 0, height = 0, nChannels = 1, nZSlices = 1, nTimepoints = 1, tileWidth = 0, tileHeight = 0;
		double pixelWidth = Double.NaN, pixelHeight = Double.NaN, zSpacing = Double.NaN, magnification = Double.NaN;
		TimeUnit timeUnit = null;
		
		// See if there is a series name embedded in the path
		int index = path.indexOf(delimiter);
		String seriesName = null;
		filePath = path;
		if (index > 0 && index < path.length()-delimiter.length() &&  !new File(path).exists()) {
			seriesName = path.substring(index+delimiter.length());
			filePath = path.substring(0, index);
		}
		
		this.path = path;
		
		// We'll need to get the metadata to parse pixel sizes etc.
	    IMetadata meta = null;

	    // Create a reader if we need to
		this.reader = cachedReaders.get(filePath);
		if (this.reader == null) {
			// Create OME-XML metadata store
		    ServiceFactory factory = new ServiceFactory();
		    OMEXMLService service = factory.getInstance(OMEXMLService.class);
		    meta = service.createOMEXMLMetadata();
			this.reader = createReader(null, filePath, meta);
		    cachedReaders.put(filePath, this.reader);
		} else {
			this.reader.setId(filePath);	
			meta = (IMetadata)this.reader.getMetadataStore();
			this.reader.setSeries(0);
		}
		
		synchronized(reader) {
		
			// Populate the image server list if we have more than one image
			int seriesIndex = -1;
			
			if (this.reader.getSeriesCount() > 1) {
				imageMap = new LinkedHashMap<>(reader.getSeriesCount());
				associatedImageMap = new LinkedHashMap<>(reader.getSeriesCount());
				if (reader.getSeriesCount() != meta.getImageCount())
					logger.error("Bio-Formats series and image counts do not match");
				for (int s = 0; s < meta.getImageCount(); s++) {
					reader.setSeries(s);
					String name = meta.getImageName(s);
					if (reader.isThumbnailSeries() || extraImageNames.contains(name.toLowerCase().trim())) {
						name = name + " (thumbnail)";
						associatedImageMap.put(name, s);
					}
					else {
						imageMap.put(name, s);
					}
					// Set this to be the series, if necessary
					if (seriesName != null && seriesName.equals(name)) {
						seriesIndex = s;
					}
					logger.debug("Adding " + meta.getImageName(s));
				}
				// If we have just one image in the image list, then reset to none - we can't switch
				if (imageMap.size() == 1 && seriesIndex < 0) {
					seriesIndex = imageMap.values().iterator().next();
					imageMap.clear();
				} else if (imageMap.size() > 1) {
					// Set default series index, if we need to
					if (seriesIndex < 0) {
						seriesIndex = imageMap.values().iterator().next();
					}
					// If we have more than one image, ensure that we have the image name correctly encoded in the path
					this.path = getSubImagePath(meta.getImageName(seriesIndex));
				}
			} else {
				if (seriesIndex < 0)
					seriesIndex = 0;
				imageMap = Collections.emptyMap();
			}
			
			
			if (seriesIndex < 0)
				throw new RuntimeException("Unable to find any non-thumbnail images within " + path);
			
			// Store the series we are actually using
			this.series = seriesIndex;
			reader.setSeries(series);
			
		    // Try getting the magnification
			if (!reader.isThumbnailSeries() && meta.getInstrumentCount() > series) {
			    try {
		    		magnification = meta.getObjectiveNominalMagnification(series, 0);
		    		if (meta.getObjectiveCount(series) > 1)
		    			logger.warn("Objective instrument count is {} - I'm not sure how to interpret this when it is != 1, check the magnification value for reasonableness", meta.getObjectiveCount(series));
			    } catch (Exception e) {
			    		logger.warn("Unable to parse magnification");
			    }
			}
			// At the time of writing, Bio-Formats does not parse the magnification from NDPI files
			// See http://openslide.org/formats/hamamatsu/
			if (Double.isNaN(magnification) && reader.getReader() instanceof NDPIReader) {
				try {
					NDPIReader ndpiReader = (NDPIReader)reader.getReader();
					IFD ifd = ndpiReader.getIFDs().get(0);
					Object o = ifd.getIFDValue(65421);
					if (o instanceof Number)
						magnification = ((Number)o).doubleValue();
				} catch (Exception e) {
					logger.warn("Problem trying to parse magnification from NDPI file", e);
				}
			}
		    		
			// Get the dimensions for the series
			int nResolutions = reader.getResolutionCount();
			int[] resWidths = new int[nResolutions];
			int[] resHeights = new int[nResolutions];
			
			for (int count = 0; count < meta.getImageCount(); count++) {
				logger.debug(meta.getImageName(count) + ": Physical pixels x: " + meta.getPixelsPhysicalSizeX(count));
			}
			
	
			Integer[] colorArray = null;
			int resLevel0 = -1;
			for (int i = 0; i < nResolutions; i++) {
				
				reader.setResolution(i);
	//			logger.debug("Series: " + series + ", Core Index: " + coreIndex);
				
				int w = reader.getSizeX();
				int h = reader.getSizeY();
				// Store the width, height & thumbnail channels if it is potentially the main image
				if (w > width && h > height) {
					width = w;
					height = h;
					resLevel0 = i;
					tileWidth = reader.getOptimalTileWidth();
					tileHeight = reader.getOptimalTileHeight();
					nChannels = reader.getSizeC();
					if (colorArray == null)
						colorArray = new Integer[nChannels];
					else if (nChannels > colorArray.length)
						colorArray = Arrays.copyOf(colorArray, nChannels);
					nZSlices = reader.getSizeZ();
					nTimepoints = reader.getSizeT();
					bpp = reader.getBitsPerPixel();
					isRGB = reader.isRGB() && bpp == 8 && nChannels == 3;
					
					// Try to read the default display colors for each channel from the file
					if (!isRGB) {
						for (int c = 0; c < nChannels; c++) {
							ome.xml.model.primitives.Color color = null;
							try {
								color = meta.getChannelColor(series, c);
							} catch (Exception e) {
								logger.warn("Unable to parse color", e);
							}
							if (color != null)
								colorArray[c] = ColorTools.makeRGBA(color.getRed(), color.getGreen(), color.getBlue(), color.getAlpha());
						}
					}
					
					// Try parsing pixel sizes in micrometers
				    try {
					    	Unit<Length> micrometers = UNITS.MICROMETER;
					    	Length xSize = meta.getPixelsPhysicalSizeX(series);
					    	Length ySize = meta.getPixelsPhysicalSizeY(series);
					    	if (xSize != null && ySize != null) {
					    		pixelWidth = xSize.value(micrometers).doubleValue();
					    		pixelHeight = ySize.value(micrometers).doubleValue();
					    	} else {
					    		pixelWidth = Double.NaN;
					    		pixelHeight = Double.NaN;			    		
					    	}
						if (nZSlices > 1) {
						    	Length z = meta.getPixelsPhysicalSizeZ(i);
						    	if (z != null)
						    		zSpacing = z.value(micrometers).doubleValue();
						    	else
						    		zSpacing = Double.NaN;
						    }
						    // TODO: Check the Bioformats TimeStamps
						    if (nTimepoints > 1) {
							    logger.warn("Time stamps read from Bioformats have not been fully verified & should not be relied upon");
							    // Here, we don't try to separate timings by z-slice & channel...
							    int lastTimepoint = -1;
							    int count = 0;
							    timePoints = new double[nTimepoints];
							    logger.debug("PLANE COUNT: " + meta.getPlaneCount(series));
							    for (int plane = 0; plane < meta.getPlaneCount(series); plane++) {
							    	int timePoint = meta.getPlaneTheT(series, plane).getValue();
							    	logger.debug("Checking " + timePoint);
							    	if (timePoint != lastTimepoint) {
							    		timePoints[count] = meta.getPlaneDeltaT(series, plane).value().doubleValue();
							    		logger.debug(String.format("Timepoint %d: %.3f seconds", count, timePoints[count]));
							    		lastTimepoint = timePoint;
							    		count++;
							    	}
							    }
							    timeUnit = TimeUnit.SECONDS;
						    }
				    } catch (Exception e) {
					    	logger.error("Error parsing metadata", e);
					    	pixelWidth = Double.NaN;
					    	pixelHeight = Double.NaN;
					    	zSpacing = Double.NaN;
					    	timePoints = null;
					    	timeUnit = null;
				    }
	
				}
				resWidths[i] = w;
				resHeights[i] = h;
			}
			
			if (resLevel0 < 0)
				throw new FormatException("Unable to find valid full-resolution image within " + path);
			
			// Loop through the series again & determine downsamples
			downsamples = new double[nResolutions];
			for (int i = 0; i < nResolutions; i++) {
				int w = resWidths[i];
				int h = resHeights[i];
				double downsampleX = (double)width / w;
				double downsampleY = (double)height / h;
				reader.setResolution(i);
				int sizeC = reader.getSizeC();
				double log2 = Math.log(2.0);
				boolean showWarning = true;
				if (GeneralTools.almostTheSame(downsampleX, downsampleY, 0.001) && sizeC == nChannels) {
					// If our downsample values are very close, average them
					downsamples[i] = (downsampleX + downsampleY) / 2;
					showWarning = false;
				}
	//				downsamples[i] = Math.round((downsampleX + downsampleY) / 2);
				else if (Math.round(downsampleX) == Math.round(downsampleY) && Math.abs(downsampleX - Math.round(downsampleX)) < 0.1 && sizeC == nChannels) {
					// If our downsample values are both close to (the same) integer, use that
					downsamples[i] = Math.round(downsampleX);
				} else {
					// Check if downsample was aiming to be close to a power of 2...
					// (This helps particularly with VSI... the true value is likely a power of 2?)
					double downsampleXlog2 = Math.log(downsampleX)/log2;
					double downsampleYlog2 = Math.log(downsampleY)/log2;
					if (Math.round(downsampleXlog2) == Math.round(downsampleYlog2)
							&& Math.abs(downsampleXlog2 - Math.round(downsampleXlog2)) < 0.1
							&& Math.abs(downsampleYlog2 - Math.round(downsampleYlog2)) < 0.1
							&& sizeC == nChannels) {
					} else {
						downsamples[i] = Double.NaN;
					}
				}
				if (showWarning)
					logger.warn("Calculated downsample values for series {} differ at resolution {}: x={} and y={} - will use value {}", series, i, downsampleX, downsampleY, downsamples[i]);
			}
			
			// Check the tile size is potentially ok
			if (tileWidth < 32 && tileWidth > 4096)
				tileWidth = -1;		
			if (tileHeight < 32 && tileHeight > 4096)
				tileHeight = -1;
			
			// Set metadata
			ImageServerMetadata.Builder builder = new ImageServerMetadata.Builder(this.path, width, height).
					setSizeC(nChannels).
					setSizeZ(nZSlices).
					setSizeT(nTimepoints).
					setPixelSizeMicrons(pixelWidth, pixelHeight).
					setZSpacingMicrons(zSpacing > 1 ? zSpacing : Double.NaN).
					setMagnification(magnification).
					setTimeUnit(timeUnit);
			if (tileWidth > 0 && tileHeight > 0)
				builder.setPreferredTileSize(tileWidth, tileHeight);
			originalMetadata = builder.build();
			
			// Set channel colors - need to wait until after metadata has been set to avoid NPE
			channelColors.clear();
			if (!isRGB) {
				for (int c = 0; c < nChannels; c++) {
					if (colorArray[c] == null)
						channelColors.add(getExtendedDefaultChannelColor(c));
					else
						channelColors.add(colorArray[c]);
				}
			}
			
		}
		
		// Bioformats can use ImageIO for JPEG decoding, and permitting the disk-based cache can slow it down... so here we turn it off
		// TODO: Document - or improve - the setting of ImageIO disk cache
		ImageIO.setUseCache(false);
		
		long endTime = System.currentTimeMillis();
		logger.debug(String.format("Initialization time: %d ms", endTime-startTime));
	}

	/**
	 * Create a new BufferedImageReader, with memoization if necessary.
	 * 
	 * @param id File path for the image.
	 * @param store Optional MetadataStore; this will be set in the reader if needed.
	 * @return the BufferedImageReader
	 * @throws FormatException
	 * @throws IOException
	 */
	private BufferedImageReader createReader(final String id, final MetadataStore store) throws FormatException, IOException {
		return createReader(null, id, store);
	}
	
	/**
	 * Create a new BufferedImageReader, with memoization if necessary.
	 * 
	 * @param cls Optionally specify a IFormatReader class if it is already known, to avoid a search.
	 * @param id File path for the image.
	 * @param store Optional MetadataStore; this will be set in the reader if needed.
	 * @return the BufferedImageReader
	 * @throws FormatException
	 * @throws IOException
	 */
	private BufferedImageReader createReader(final Class<? extends IFormatReader> cls, final String id, final MetadataStore store) throws FormatException, IOException {
		IFormatReader imageReader;
		if (cls != null) {
			ClassList<IFormatReader> list = new ClassList<>(IFormatReader.class);
			list.addClass(reader.getReader().getClass());
			imageReader = new ImageReader(list);
		} else
			imageReader = new ImageReader();
		
		imageReader.setFlattenedResolutions(false);
		
		Memoizer memoizer = null;
		if (options.requestMemoization()) {
			String pathMemoization = options.getPathMemoization();
			long memoizationMinTimeMillis = 100L;
			if (pathMemoization != null && !pathMemoization.trim().isEmpty()) {
				File dir = new File(pathMemoization);
				if (dir.isDirectory())
					memoizer = new Memoizer(imageReader, memoizationMinTimeMillis, dir);
				else {
					logger.warn("Memoization directory '{}' not found - will default to image directory", pathMemoization);
					memoizer = new Memoizer(imageReader, memoizationMinTimeMillis);
				}
			} else
				memoizer = new Memoizer(imageReader, memoizationMinTimeMillis);
			imageReader = memoizer;
		}
		
		if (store != null)
			imageReader.setMetadataStore(store);
		
		if (id != null) {
			imageReader.setId(id);
			if (memoizer != null) {
				File fileMemo = ((Memoizer)imageReader).getMemoFile(filePath);
				memoizationFileSize = fileMemo == null ? 0L : fileMemo.length();
				if (memoizationFileSize == 0L)
					logger.debug("No memoization file generated for {}", id);
				else
					logger.debug(String.format("Generating memoization file %s (%.2f MB)", fileMemo.getAbsolutePath(), memoizationFileSize/1024.0/1024.0));			
			}
		}
		return BufferedImageReader.makeBufferedImageReader(imageReader);
	}
	
		
	/**
	 * Returns true if the reader accepts parallel tile requests, without synchronization.
	 * 
	 * This is true if parallelization is requested, and any memoization file is less than 10 MB.
	 * The idea is that larger memoization files indicate more heavyweight readers, and these need 
	 * to be kept restricted.
	 * 
	 * @return
	 */
	public boolean willParallelize() {
		return options.requestParallelization() && getWidth() > 8192 && getHeight() > 8192 && memoizationFileSize < 1024L*1024L * 10L;
	}
	
	/**
	 * Get an IFormatReader for the current thread, storing the result in a map if needed.
	 * 
	 * @return
	 */
	private IFormatReader getCurrentThreadImageReader() {
		// Here, to avoid requiring synchronization and to avoid trouble with the reader state,
		// we create readers for each thread
		// This doesn't appear to be too expensive if memoization is turned on, and the readers are sufficiently lightweight.
		// See https://docs.openmicroscopy.org/bio-formats/5.7.2/developers/matlab-dev.html#improving-reading-performance
		IFormatReader ifReader = map.get(Thread.currentThread());
		if (ifReader == null) {
			try {
				ifReader = createReader(reader.getReader().getClass(), filePath, null);
				map.put(Thread.currentThread(), ifReader);
			} catch (ClosedByInterruptException e) {
				logger.error(e.getLocalizedMessage());
			} catch (Exception e) {
				logger.error(e.getLocalizedMessage());
			}
		}
		try {
			// Make extra sure the ID is set
			if (ifReader == null) {
				ifReader = createReader(filePath, null);
			}
			if (ifReader.getCurrentFile() == null)
				ifReader.setId(filePath);			
		} catch (Exception e) {
			logger.error("Problem getting an image reader", e);
		}
		return ifReader;
	}
	
	/**
	 * Get a BufferedImageReader for use by the current thread.
	 * 
	 * If willParallelize() returns false, then the global reader will be provided.
	 * 
	 * @return
	 */
	private BufferedImageReader getBufferedImageReader() {
		IFormatReader ifReader = willParallelize() ? getCurrentThreadImageReader() : reader;
		return BufferedImageReader.makeBufferedImageReader(ifReader);
	}
	
	
	@Override
	public PathImage<BufferedImage> readRegion(RegionRequest request) {
		return new PathBufferedImage(this, request, readBufferedImage(request));
	}

	@Override
	public BufferedImage readBufferedImage(RegionRequest request) {
		
		int resolution = ServerTools.getClosestDownsampleIndex(getPreferredDownsamples(), request.getDownsample());
		double downsampleFactor = request.getDownsample();
		double downsampleForSeries = getPreferredDownsamples()[resolution];
				
		/*
		 * For downsampled levels, Aperio .svs files appear to divide by 4, 32 etc. and floor (i.e. not round).
		 * Consequently their downsample ratios are not round numbers.
		 * The numbers reported in the metadata appear to be the mean of the width & height ratios for each level
		 * (i.e. full_res_width / level width and full_res_height / level_height).
		 * This can lead to values like 4.000112795792701... close to 4, but not quite there.
		 * 
		 * This is... troublesome, because it seems Bioformats BufferedImageReader can put some black borders
		 * on image tiles when their coordinates are slightly off from the real TIFF tiles.
		 * 
		 * However by using ceil (rather than round) in the scaling below it seems to be ok.
		 * Still, it is something to be wary of, and it may need revisiting - especially with other formats.
		 * 
		 */
		
		// Adjust coordinates if we are downsampling
		Rectangle region2;
		if (downsampleForSeries > 1) {
			region2 = new Rectangle(
					(int)Math.ceil(request.getX() / downsampleForSeries),
					(int)Math.ceil(request.getY() / downsampleForSeries),
					(int)Math.ceil(request.getWidth() / downsampleForSeries),
					(int)Math.ceil(request.getHeight() / downsampleForSeries));
		} else {
			region2 = AwtTools.getBounds(request);
		}
		
		BufferedImage img;
		
		BufferedImageReader ipReader = getBufferedImageReader();
		synchronized(ipReader) {
			try {
				if (ipReader.getCurrentFile() == null)
					ipReader.setId(filePath);
				ipReader.setSeries(series);
				ipReader.setResolution(resolution);
				
				// Ensure the region coordinates are within range
				region2.x = region2.x < 0 ? 0 : region2.x;
				region2.y = region2.y < 0 ? 0 : region2.y;
				region2.width = region2.x + region2.width >= ipReader.getSizeX() ? ipReader.getSizeX() - region2.x : region2.width;
				region2.height = region2.y + region2.height >= ipReader.getSizeY() ? ipReader.getSizeY() - region2.y : region2.height;
				
				// Determine the final required size - which may or may not be the same
				int finalWidth, finalHeight;
				boolean resizeRequired;
				if (GeneralTools.almostTheSame(downsampleForSeries, downsampleFactor, 0.001)) {
					finalWidth = request.getWidth();
					finalHeight = request.getHeight();
					resizeRequired = false;
				} else {
					finalWidth = (int)(request.getWidth() / downsampleFactor + .5);
					finalHeight = (int)(request.getHeight() / downsampleFactor + .5);
					resizeRequired = true;
				}

				// Read the image - or at least the first channel
				int ind = ipReader.getIndex(request.getZ(), 0, request.getT());
				img = ipReader.openImage(ind, region2.x, region2.y, region2.width, region2.height);
				
				// Resize if we need to
				if (resizeRequired)
					img = resize(img, finalWidth, finalHeight, isRGB);
				
				// Single-channel & RGB images are straightforward... nothing more to do
				if (ipReader.isRGB() || nChannels() == 1)
					return img;
				
				// If we have multiple channels, merge them
				BufferedImage[] images = new BufferedImage[nChannels()];
				images[0] = img;
				for (int c = 1; c < nChannels(); c++) {
					ind = ipReader.getIndex(request.getZ(), c, request.getT());
					img = ipReader.openImage(ind, region2.x, region2.y, region2.width, region2.height);
					if (resizeRequired)
						img = resize(img, finalWidth, finalHeight, isRGB);
					images[c] = img;
				}
				return AWTImageTools.mergeChannels(images);

			} catch (Exception e) {
				logger.error("Error reading image region " + request + " for image size " + ipReader.getSizeX() + " x " + ipReader.getSizeY(), e);
			}
		}
		return null;
	}
	
	/**
	 * Resize the image to have the requested width/height.
	 * 
	 * @param img
	 * @param finalWidth
	 * @param finalHeight
	 * @return
	 */
	private BufferedImage resize(final BufferedImage img, final int finalWidth, final int finalHeight, final boolean isRGB) {
		// RGB can generally be converted more easily
		if (isRGB) {
			try {
				BufferedImage img2 = new BufferedImage(finalWidth, finalHeight, img.getType());
				Graphics2D g2d = img2.createGraphics();
				g2d.drawImage(img, 0, 0, finalWidth, finalHeight, null);
				g2d.dispose();
				return img2;
			} catch (Exception e) {
				logger.debug("Error rescaling (supposedly) RGB image {}, will default to slower rescaling: {}", img, e.getLocalizedMessage());
			}
		}
		
		// Get the pixels
		float[] pixels = img.getRaster().getPixels(0, 0, img.getWidth(), img.getHeight(), (float[])null);
		double xScale = (double)img.getWidth() / finalWidth;
		double yScale = (double)img.getHeight() / finalHeight;
		
		// Perform rescaling with nearest neighbor interpolation
		// TODO: Consider 'better' forms of interpolation
		float[] pixelsNew = new float[finalWidth*finalHeight];
		int w = img.getWidth();
		int h = img.getHeight();
		for (int y = 0; y < finalHeight; y++) {
			int row = (int)(y * yScale + 0.5);
			if (row >= h)
				row = h-1;
			for (int x = 0; x < finalWidth; x++) {
				int col = (int)(x * xScale + 0.5);
				if (col >= w)
					col = w-1;
				int ind = row*img.getWidth() + col;
				pixelsNew[y*finalWidth + x] = pixels[ind];
			}			
		}
		
		// Create an image with the same ColorModel / data type as the original
		WritableRaster raster = img.getColorModel().createCompatibleWritableRaster(finalWidth, finalHeight);
		raster.setSamples(0, 0, finalWidth, finalHeight, 0, pixelsNew);
		return new BufferedImage(img.getColorModel(), raster, img.getColorModel().isAlphaPremultiplied(), null);
		
//		// Warning!  This doesn't actually work!  It performs some unwelcome rescaling of pixel intensities
//		logger.warn("Resizing not implemented properly for images with type {} - pixel values will be surreptitiously rescaled", img.getType());
//		return AWTImageTools.scale(img, finalWidth, finalHeight, false);
	}
	
	
	@Override
	public String getServerType() {
		return "Bio-Formats";
	}
	
	@Override
	public boolean isRGB() {
		return isRGB;
	}
	
	@Override
	public double[] getPreferredDownsamples() {
		return downsamples.clone();
	}

	@Override
	public synchronized void close() {
		super.close();
//		if (reader != null) {
//			try {
//				reader.close();
//			} catch (IOException e) {
//				e.printStackTrace();
//			}
//		}

		for (IFormatReader reader : map.values()) {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	@Override
	public BufferedImage getBufferedThumbnail(int maxWidth, int maxHeight, int zPosition) {
		BufferedImage img = super.getBufferedThumbnail(maxWidth, maxHeight, zPosition);
		if (isRGB())
			return img;
		return AWTImageTools.autoscale(img);
	}
	
	@Override
	public double getTimePoint(int ind) {
		if (nTimepoints() == 0)
			return 0;
		return timePoints[ind];
	}

	@Override
	public List<String> getSubImageList() {
		if (imageMap == null || imageMap.isEmpty())
			return Collections.emptyList();
		return new ArrayList<>(imageMap.keySet());
	}

	@Override
	public String getDisplayedImageName() {
		return getShortServerName();
	}

	@Override
	public boolean containsSubImages() {
		return imageMap != null && !imageMap.isEmpty();
	}


	@Override
	public boolean usesBaseServer(ImageServer<?> server) {
		return this == server;
	}


	@Override
	public int getBitsPerPixel() {
		return bpp;
	}
	
	
	@Override
	public Integer getDefaultChannelColor(int channel) {
		if (isRGB)
			return getDefaultRGBChannelColors(channel);
		Integer color = null;
		if (channel >= 0 && channel <= channelColors.size())
			 color = channelColors.get(channel);
		if (color == null)
			color = getExtendedDefaultChannelColor(channel);
		return color;
	}
	
	@Override
	public List<String> getAssociatedImageList() {
		if (associatedImageMap == null || associatedImageMap.isEmpty())
			return Collections.emptyList();
		return new ArrayList<>(associatedImageMap.keySet());
	}

	@Override
	public BufferedImage getAssociatedImage(String name) {
		if (associatedImageMap == null || !associatedImageMap.containsKey(name))
			throw new IllegalArgumentException("No associated image with name '" + name + "' for " + getPath());
		BufferedImageReader reader = getBufferedImageReader();
		synchronized (reader) {
			int series = reader.getSeries();
			try {
				reader.setSeries(associatedImageMap.get(name));
				if (reader.getResolutionCount() > 0) {
					reader.setResolution(0);
				}
				// TODO: Handle color transforms here, or in the display of labels/macro images - in case this isn't RGB
				BufferedImage img = reader.openImage(reader.getIndex(0, 0, 0), 0, 0, reader.getSizeX(), reader.getSizeY());
				return img;
			} catch (Exception e) {
				logger.error("Error reading associated image" + name, e);
			} finally {
				reader.setSeries(series);
			}
		}
		return null;
	}
	
	
	@Override
	public File getFile() {
		File file = new File(filePath);
		if (file.exists())
			return file;
		return null;
	}


	@Override
	public String getSubImagePath(String imageName) {
		if (imageMap.containsKey(imageName))
			return filePath + delimiter + imageName;
		throw new IllegalArgumentException(toString() + " does not contain sub-image with name " + imageName);
	}


	@Override
	public ImageServerMetadata getMetadata() {
		return userMetadata == null ? originalMetadata : userMetadata;
	}


	@Override
	public ImageServerMetadata getOriginalMetadata() {
		return originalMetadata;
	}

	@Override
	public void setMetadata(ImageServerMetadata metadata) {
		if (!originalMetadata.isCompatibleMetadata(metadata))
			throw new RuntimeException("Specified metadata is incompatible with original metadata for " + this);
		userMetadata = metadata;
	}
	
}