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

import java.awt.image.BufferedImage;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import qupath.lib.images.servers.FileFormatInfo.ImageCheckType;

/**
 * Builder for ImageServers that make use of the Bio-Formats library.
 * 
 * @author Pete Bankhead
 *
 */
public class BioFormatsServerBuilder implements ImageServerBuilder<BufferedImage> {
	
	final private static Logger logger = LoggerFactory.getLogger(BioFormatsServerBuilder.class);

	@Override
	public ImageServer<BufferedImage> buildServer(String path) {
		try {
			return new BioformatsImageServer(path);
		} catch (Exception e) {
			logger.error("Unable to open {}", path, e);
		}
		return null;
	}

	@Override
	public float supportLevel(String path, ImageCheckType type, Class<?> cls) {
		if (cls != BufferedImage.class)
			return 0;
		switch (type) {
		case TIFF_2D_RGB:
			return 3;
		case TIFF_IMAGEJ:
			return 3;
		case TIFF_OTHER:
			return 2;
		case UNKNOWN:
			return 2;
		case URL:
			return 0;
		default:
			return 2;
		}
	}

	@Override
	public String getName() {
		return "Bio-Formats";
	}

	@Override
	public String getDescription() {
		return "Image server using the Bio-Formats library";
	}
	
}