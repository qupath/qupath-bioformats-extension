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
import java.util.concurrent.TimeUnit;

import javax.imageio.ImageIO;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ome.units.quantity.Length;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.ClassList;
import loci.formats.FormatException;
import loci.formats.IFormatReader;
import loci.formats.ImageReader;
import loci.formats.gui.AWTImageTools;
import loci.formats.gui.BufferedImageReader;
import loci.formats.in.NDPIReader;
import loci.formats.meta.IMetadata;
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
public class BioformatsImageServer extends AbstractImageServer<BufferedImage> {
	
	private static final Logger logger = LoggerFactory.getLogger(BioformatsImageServer.class);
	
	private ImageServerMetadata originalMetadata;
	private ImageServerMetadata userMetadata;

	private String path;

	private String filePath; // Path to the base image file - will be the same as path, unless the path encodes the name of a specific series, in which case this refers to the file without the series included
	
	private ImageReader reader;
	private int resLevel0 = -1;
	private double[] downsamples;
	private boolean isRGB = true;
	private int bpp = 0;
	
	// For image lists, store name & series number
	private Map<String, Integer> imageMap = null;
	private Map<String, Integer> associatedImageMap = null;

	private List<Integer> channelColors = new ArrayList<>();

	private double[] timePoints = null;
	
	private String delimiter = "::"; // Delimiter between the file path and any sub-image names

	private IMetadata meta;
	private int series = 0; // This is effectively the image (there might be more than one in the file)
		
	public int getCurrentSeries() {
		return series;
	}
	

	public BioformatsImageServer(final String path) throws FormatException, IOException, DependencyException, ServiceException {
		this(new ImageReader(), path);
	}


	BioformatsImageServer(final ImageReader reader, final String path) throws FormatException, IOException, DependencyException, ServiceException {

		super();

		// Warning, because my interpretation of what Bio-Formats is doing may be wrong...
		// leading to some unfortunate bugs
		logger.warn("QuPath Bio-Formats extension is in beta!  Be watchful for bugs...");

		this.path = path;
		this.reader = reader;
		
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

		// Zip files tend to be slow for Bioformats to parse... so best not try
		if (path.endsWith(".zip"))
			throw new FormatException("ZIP not supported");
		
		// We want each series to be a different image - not to have resolutions 'flattened' into series
	    reader.setFlattenedResolutions(false);
		
		// Create OME-XML metadata store (so we can get pixel sizes)
	    ServiceFactory factory = new ServiceFactory();
	    OMEXMLService service = factory.getInstance(OMEXMLService.class);
	    meta = service.createOMEXMLMetadata();
	    reader.setMetadataStore(meta);
//	    reader.setOriginalMetadataPopulated(true);

	    // Set the path ID
		reader.setId(filePath);
		
		// Populate the image server list if we have more than one image
		int seriesIndex = -1;
		if (reader.getSeriesCount() > 1) {
		    logger.debug("WARNING! Metadata parsed from files containing multiple images using Bioformats may be incorrect!  Please check pixel sizes.");

			imageMap = new LinkedHashMap<>(reader.getSeriesCount());
			associatedImageMap = new LinkedHashMap<>(reader.getSeriesCount());
			logger.debug("Reader series count: " + reader.getSeriesCount() + ", Reader image count: " + reader.getImageCount() + ", Metadata image count: " + meta.getImageCount());
			if (reader.getSeriesCount() != meta.getImageCount())
				logger.error("WARNING! Bio-Formats series and image counts do not match");
			for (int s = 0; s < meta.getImageCount(); s++) {
				reader.setSeries(s);
				String name;
				if (reader.isThumbnailSeries()) {
					name = meta.getImageName(s) + " (thumbnail)";
					associatedImageMap.put(name, s);
				}
				else {
					name = meta.getImageName(s);
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
				// Set default series index
				// This code is rather terrible, but it strives to accept the first image not called 'overview' (which appears to be the first image in some VSI files)
				if (seriesIndex < 0) {
					for (String key : imageMap.keySet()) {
						if (!key.toLowerCase().equals("overview")) {
							seriesIndex = imageMap.get(key);
							break;
						}
					}
					if (seriesIndex < 0)
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
	    			logger.debug("WARNING! Objective instrument count is " + meta.getObjectiveCount(series) + " - I'm not sure how to interpret this when it is != 1");
		    } catch (Exception e) {
		    	logger.debug("Unable to parse magnification");
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
				e.printStackTrace();
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
		for (int i = 0; i < nResolutions; i++) {
			
			reader.setResolution(i);
			int imageIndex = series; // TODO: Fix this!  It appears to be wrong sometimes...?
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
				isRGB = reader.isRGB();
				bpp = reader.getBitsPerPixel();
				
				// Try to read the default display colors for each channel from the file
				if (!isRGB) {
					for (int c = 0; c < nChannels; c++) {
						ome.xml.model.primitives.Color color = meta.getChannelColor(imageIndex, c);
						if (color != null)
							colorArray[c] = ColorTools.makeRGBA(color.getRed(), color.getGreen(), color.getBlue(), color.getAlpha());
//						if (color != null)
//							channelColors.add(ColorTools.makeRGBA(color.getRed(), color.getGreen(), color.getBlue(), color.getAlpha()));
//						else
//							channelColors.add(getExtendedDefaultChannelColor(c));
					}
				}
				
				// Try parsing pixel sizes
			    try {
			    	// TODO: Check if the assumption that this values are microns actually holds in all cases...
			    	Length xSize = meta.getPixelsPhysicalSizeX(imageIndex);
			    	Length ySize = meta.getPixelsPhysicalSizeY(imageIndex);
			    	if (xSize != null && ySize != null) {
			    		pixelWidth = xSize.value().doubleValue();
			    		pixelHeight = ySize.value().doubleValue();
			    	} else {
			    		pixelWidth = Double.NaN;
			    		pixelHeight = Double.NaN;			    		
			    	}
				    if (nZSlices > 1) {
				    	Length z = meta.getPixelsPhysicalSizeZ(i);
				    	if (z != null)
				    		zSpacing = z.value().doubleValue();
				    	else
				    		zSpacing = Double.NaN;
				    }
				    // TODO: Check the Bioformats TimeStamps
				    if (nTimepoints > 1) {
					    logger.debug("WARNING! Time stamps read from Bioformats have not been fully verified & may be incorrect");
					    // Here, we don't try to separate timings by z-slice & channel...
					    int lastTimepoint = -1;
					    int count = 0;
					    timePoints = new double[nTimepoints];
					    logger.debug("PLANE COUNT: " + meta.getPlaneCount(imageIndex));
					    for (int plane = 0; plane < meta.getPlaneCount(imageIndex); plane++) {
					    	int timePoint = meta.getPlaneTheT(imageIndex, plane).getValue();
					    	logger.debug("Checking " + timePoint);
					    	if (timePoint != lastTimepoint) {
					    		timePoints[count] = meta.getPlaneDeltaT(imageIndex, plane).value().doubleValue();
					    		logger.debug(String.format("Timepoint %d: %.3f seconds", count, timePoints[count]));
					    		lastTimepoint = timePoint;
					    		count++;
					    	}
					    }
					    timeUnit = TimeUnit.SECONDS;
				    }
			    } catch (Exception e) {
			    	e.printStackTrace();
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
			if (GeneralTools.almostTheSame(downsampleX, downsampleY, 0.001) && sizeC == nChannels)
				downsamples[i] = (downsampleX + downsampleY) / 2;
//				downsamples[i] = Math.round((downsampleX + downsampleY) / 2);
			else
				downsamples[i] = Double.NaN;
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
		
		// Bioformats can use ImageIO for JPEG decoding, and permitting the disk-based cache can slow it down... so here we turn it off
		// TODO: Document - or improve - the setting of ImageIO disk cache
		ImageIO.setUseCache(false);
	}
	
	
	
	private Map<Thread, IFormatReader> map = new HashMap<>();
	
	
	public IFormatReader getImageReader() {
		// I admit... I'm a bit vague on how Bioformats works
		// Here, to avoid requiring synchronization and to avoid trouble with the reader state,
		// we create readers for each thread... it may be this is a Very Bad Idea, and it certainly
		// seems to require too much metadata parsing - but it's the best I have so far...
		IFormatReader ifReader = map.get(Thread.currentThread());
		if (ifReader == null) {
//			synchronized(this) {
//				try {
//					ifReader = reader.getReader().getClass().newInstance();				
//					ifReader.setFlattenedResolutions(false);
//					ifReader.setId(path);
//					map.put(Thread.currentThread(), ifReader);
//					return ifReader;
//				} catch (Exception e) {
//					e.printStackTrace();
//				}
//			}
			try {
				ClassList<IFormatReader> list = new ClassList<>(IFormatReader.class);
				list.addClass(reader.getReader().getClass());
				ifReader = new ImageReader(list);
				ifReader.setFlattenedResolutions(false);
				ifReader.setId(filePath);
				map.put(Thread.currentThread(), ifReader);
			} catch (ClosedByInterruptException e) {
				logger.error(e.getLocalizedMessage());
			} catch (Exception e) {
				logger.error(e.getLocalizedMessage());
//				e.printStackTrace();
			}
		}
		try {
			// Make sure the ID is set
			if (ifReader == null) {
				ifReader = new ImageReader();
				ifReader.setFlattenedResolutions(false);
				ifReader.setId(filePath);
			}
			if (ifReader.getCurrentFile() == null)
				ifReader.setId(filePath);			
		} catch (Exception e) {
			e.printStackTrace();
		}
		return ifReader;
	}
	
	
	private BufferedImageReader getBufferedImageReader() {
		IFormatReader ifReader = getImageReader();
		return BufferedImageReader.makeBufferedImageReader(ifReader);
	}
	
	
	@Override
	public PathImage<BufferedImage> readRegion(RegionRequest request) {
		return new PathBufferedImage(this, request, readBufferedImage(request));
	}

	@Override
	public BufferedImage readBufferedImage(RegionRequest request) {
		
		BufferedImageReader ipReader = getBufferedImageReader();
		
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
		
//		synchronized(this) 
		{
			try {
				if (ipReader.getCurrentFile() == null)
					ipReader.setId(filePath);
				ipReader.setSeries(series);
				ipReader.setResolution(resolution);
				
				// Ensure the region coordinates are within range
				region2 = region2.intersection(new Rectangle(0, 0, ipReader.getSizeX(), ipReader.getSizeY()));
				
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
				logger.error("Error reading image region", e);
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
		if (isRGB || img.getType() == BufferedImage.TYPE_BYTE_GRAY) {
			BufferedImage img2 = new BufferedImage(finalWidth, finalHeight, BufferedImage.TYPE_INT_RGB);
			Graphics2D g2d = img2.createGraphics();
			g2d.drawImage(img, 0, 0, finalWidth, finalHeight, null);
			g2d.dispose();
			return img2;
		}
		
		// Get the pixels
		float[] pixels = AWTImageTools.getFloats(img)[0];
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
		if (reader != null) {
			try {
				reader.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
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
//		ImageServer server = getSubImageServer(name);
		ImageServer<BufferedImage> server;
		try {
			server = new BioformatsImageServer(filePath + delimiter + name);
			BufferedImage img = server.getBufferedThumbnail(2048, 2048, 0);
			if (!server.usesBaseServer(this))
				server.close();
			return img;
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
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
		throw new RuntimeException(toString() + " does not contain sub-image with name " + imageName);
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
