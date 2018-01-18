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

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.StringTokenizer;

import javafx.beans.property.BooleanProperty;
import javafx.beans.property.StringProperty;
import qupath.lib.gui.QuPathGUI;
import qupath.lib.gui.extensions.QuPathExtension;
import qupath.lib.gui.panels.PreferencePanel;
import qupath.lib.gui.prefs.PathPrefs;

/**
 * A QuPath extension that adds options relating to the BioFormatsImageServer to the main QuPath preference pane.
 * 
 * @author Pete Bankhead
 */
public class BioFormatsOptionsExtension implements QuPathExtension {

	@Override
	public void installExtension(QuPathGUI qupath) {
		
		BioFormatsServerOptions options = BioFormatsServerOptions.getInstance();
		
		// Create persistent properties
		BooleanProperty enableBioformats = PathPrefs.createPersistentPreference("bfEnableBioformats", options.bioformatsEnabled());
		BooleanProperty useParallelization = PathPrefs.createPersistentPreference("bfUseParallization", options.requestParallelization());
		BooleanProperty useMemoization = PathPrefs.createPersistentPreference("bfUseMemoization", options.requestMemoization());
		
		StringProperty pathMemoization = PathPrefs.createPersistentPreference("bfPathMemoization", options.getPathMemoization());
		StringProperty useExtensions = PathPrefs.createPersistentPreference("bfUseAlwaysExtensions", String.join(" ", options.getUseAlwaysExtensions()));
		StringProperty skipExtensions = PathPrefs.createPersistentPreference("bfSkipAlwaysExtensions", String.join(" ", options.getSkipAlwaysExtensions()));
		
		// Set options using any values previously stored
		options.setPathMemoization(pathMemoization.get());
		options.setBioformatsEnabled(enableBioformats.get());
		options.setRequestParallelization(useParallelization.get());
		options.setRequestMemoization(useMemoization.get());
		fillCollectionWithTokens(useExtensions.get(), options.getUseAlwaysExtensions());
		fillCollectionWithTokens(skipExtensions.get(), options.getSkipAlwaysExtensions());

		// Listen for property changes
		enableBioformats.addListener((v, o, n) -> options.setBioformatsEnabled(n));
		useParallelization.addListener((v, o, n) -> options.setRequestParallelization(n));
		useMemoization.addListener((v, o, n) -> options.setRequestMemoization(n));
		
		pathMemoization.addListener((v, o, n) -> options.setPathMemoization(n));
		useExtensions.addListener((v, o, n) -> fillCollectionWithTokens(n, options.getUseAlwaysExtensions()));
		skipExtensions.addListener((v, o, n) -> fillCollectionWithTokens(n, options.getSkipAlwaysExtensions()));
		
		// Add preferences to QuPath GUI
		PreferencePanel prefs = QuPathGUI.getInstance().getPreferencePanel();
		prefs.addPropertyPreference(enableBioformats, Boolean.class, "Enable Bio-Formats", "Bio-Formats", "Allow QuPath to use Bio-Formats for image reading");
		prefs.addPropertyPreference(useParallelization, Boolean.class, "Enable Bio-Formats parallelization", "Bio-Formats", "Enable reading image tiles in parallel when using Bio-Formats");
		prefs.addPropertyPreference(useMemoization, Boolean.class, "Enable Bio-Formats memoization", "Bio-Formats", "Allow Bio-Formats to use memoization to improve performance when reopening large or complex images");
		
		prefs.addDirectoryPropertyPreference(pathMemoization, "Bio-Formats memoization directory", "Bio-Formats",
				"Choose directory where Bio-Formats should write cache files for memoization; by default the directory where the image is stored will be used");
		prefs.addPropertyPreference(useExtensions, String.class, "Always use Bio-Formats for image extensions", "Bio-Formats", 
				"Request that Bio-Formats is always the file reader used for images with specific extensions; enter as a list with spaces between each entry");
		prefs.addPropertyPreference(skipExtensions, String.class, "Never use Bio-Formats for image extensions", "Bio-Formats", 
				"Request that Bio-Formats is never the file reader used for images with specific extensions; enter as a list with spaces between each entry");

	}

	private static void fillCollectionWithTokens(String text, Collection<String> collection) {
		fillCollectionWithTokens(new StringTokenizer(text), collection);
	}

	private static void fillCollectionWithTokens(StringTokenizer tokenizer, Collection<String> collection) {
		List<String> list = new ArrayList<>();
		while (tokenizer.hasMoreTokens())
			list.add(tokenizer.nextToken());
		collection.clear();
		collection.addAll(list);
	}

	@Override
	public String getName() {
		return "Bio-Formats server options";
	}

	@Override
	public String getDescription() {
		return "Installs options for the Bio-Formats image server in the QuPath preference pane";
	}

}
