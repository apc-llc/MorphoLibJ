/*-
 * #%L
 * Mathematical morphology library and plugins for ImageJ/Fiji.
 * %%
 * Copyright (C) 2014 - 2017 INRA.
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Lesser Public License for more details.
 * 
 * You should have received a copy of the GNU General Lesser Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/lgpl-3.0.html>.
 * #L%
 */
package inra.ijpb.plugins;

import static org.junit.Assert.assertEquals;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ByteProcessor;
import ij.process.ImageConverter;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import inra.ijpb.watershed.Watershed;

import java.util.Map;
import java.util.HashMap;

import org.junit.Test;

public class TestWatershed2D {

	long xorshift_state = 1;

	/**
	 * Random number generator to generate random mask.
	 */
	long xorshift()
	{
		long x = xorshift_state;
		x ^= (x << 21);
		x ^= (x >>> 35);
		x ^= (x << 4);
		xorshift_state = x;
		return x;
	}
	
	ImagePlus randomImage(int width, int height)
	{
		ImageProcessor ip = new ByteProcessor( width, height );
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++)
				ip.setf(x, y, xorshift() % 2);
		}

		return new ImagePlus("", ip);
	}

	/**
	 * Test 2D watershed over the same image stored as 8, 16 and 32-bit.
	 */
	@Test
	public void testWatershed()
	{
		ImagePlus input = IJ.openImage( TestWatershed2D.class.getResource( "/files/grains.tif" ).getFile() );

		final ImagePlus copy = input.duplicate();
		
		int[] connectivityValues = new int[]{ 4, 8 };

		for( int connectivity : connectivityValues )
		{
			input = copy;

			ImageProcessor inputIP = input.getChannelProcessor();

			ImagePlus mask = randomImage( inputIP.getWidth(), inputIP.getHeight() );

			ImageProcessor maskIP = mask.getChannelProcessor();

			ImageProcessor resultIP8bit = Watershed.computeWatershed( inputIP, maskIP, connectivity );
			ImagePlus result8bit = new ImagePlus("result", resultIP8bit);

			ImagePlus reference = IJ.openImage( TestWatershed2D.class.getResource(
				"/files/grains_watershed_" + connectivity + ".tif" ).getFile() );
			
			assertEquals( "Different result compared to reference image (connectivity = " + connectivity + ")",
				0, diffImagePlus( result8bit, reference ) );

			ImageConverter ic = new ImageConverter( input );
			ic.convertToGray16();

			ImageProcessor resultIP16bit = Watershed.computeWatershed( inputIP, maskIP, connectivity );
			ImagePlus result16bit = new ImagePlus("result", resultIP16bit);

			assertEquals( "Different results for 8 and 16 bit images (connectivity = " + connectivity + ")", 
				0, diffImagePlus( result8bit, result16bit ) );

			ic = new ImageConverter( input.duplicate() );
			ic.convertToGray32();

			ImageProcessor resultIP32bit = Watershed.computeWatershed( inputIP, maskIP, connectivity );
			ImagePlus result32bit = new ImagePlus("result", resultIP32bit);

			assertEquals( "Different results for 8 and 32 bit images (connectivity = " + connectivity + ")",
				0, diffImagePlus( result8bit, result32bit ) );
		}
	}

	/**
	 * Benchmark 2D watershed over the same image stored as 8, 16 and 32-bit.
	 */
	@Test
	public void benchmarkWatershed()
	{
		Map<String, Double> means = new HashMap<>();
		means.put("[16][4][8]", 3.87890625);
		means.put("[16][8][8]", 2.375);
		means.put("[16][4][16]", 3.87890625);
		means.put("[16][8][16]", 2.375);
		means.put("[16][4][32]", 3.87890625);
		means.put("[16][8][32]", 2.375);
		means.put("[32][4][8]", 14.978515625);
		means.put("[32][8][8]", 9.8740234375);
		means.put("[32][4][16]", 14.978515625);
		means.put("[32][8][16]", 9.8740234375);
		means.put("[32][4][32]", 14.978515625);
		means.put("[32][8][32]", 9.8740234375);
		means.put("[64][4][8]", 65.816162109375);
		means.put("[64][8][8]", 36.99169921875);
		means.put("[64][4][16]", 65.816162109375);
		means.put("[64][8][16]", 36.99169921875);
		means.put("[64][4][32]", 65.816162109375);
		means.put("[64][8][32]", 36.99169921875);
		means.put("[128][4][8]", 262.84185791015625);
		means.put("[128][8][8]", 143.35784912109375);
		means.put("[128][4][16]", 262.84185791015625);
		means.put("[128][8][16]", 143.35784912109375);
		means.put("[128][4][32]", 262.84185791015625);
		means.put("[128][8][32]", 143.35784912109375);
		means.put("[256][4][8]", 1032.5812225341797);
		means.put("[256][8][8]", 565.5418853759766);
		means.put("[256][4][16]", 1032.5812225341797);
		means.put("[256][8][16]", 565.5418853759766);
		means.put("[256][4][32]", 1032.5812225341797);
		means.put("[256][8][32]", 565.5418853759766);
		means.put("[512][4][8]", 4161.341079711914);
		means.put("[512][8][8]", 2318.1643676757812);
		means.put("[512][4][16]", 4161.341079711914);
		means.put("[512][8][16]", 2318.1643676757812);
		means.put("[512][4][32]", 4161.341079711914);
		means.put("[512][8][32]", 2318.1643676757812);
		means.put("[1024][4][8]", 16581.678923606873);
		means.put("[1024][8][8]", 9277.968538284302);
		means.put("[1024][4][16]", 16581.678923606873);
		means.put("[1024][8][16]", 9277.968538284302);
		means.put("[1024][4][32]", 16581.678923606873);
		means.put("[1024][8][32]", 9277.968538284302);
		means.put("[2048][4][8]", 66414.02896332741);
		means.put("[2048][8][8]", 36933.99206995964);
		means.put("[2048][4][16]", 66414.02896332741);
		means.put("[2048][8][16]", 36933.99206995964);
		means.put("[2048][4][32]", 66414.02896332741);
		means.put("[2048][8][32]", 36933.99206995964);
		means.put("[4096][4][8]", 264885.5006724596);
		means.put("[4096][8][8]", 147389.29292684793);
		means.put("[4096][4][16]", 264885.5006724596);
		means.put("[4096][8][16]", 147389.29292684793);
		means.put("[4096][4][32]", 264885.5006724596);
		means.put("[4096][8][32]", 147389.29292684793);

		for (int dim = 16; dim <= 4096; dim *= 2)
		{
			ImagePlus input = randomImage( dim, dim );

			final ImagePlus copy = input.duplicate();
		
			int[] connectivityValues = new int[]{ 4, 8 };

			ImagePlus mask = randomImage( dim, dim );

			for( int connectivity : connectivityValues )
			{
				input = copy;

				ImageProcessor inputIP = input.getChannelProcessor();
				
				ImageProcessor maskIP = mask.getChannelProcessor();

				ImageProcessor resultIP8bit;
				{
					final long t0 = System.currentTimeMillis();
					resultIP8bit = Watershed.computeWatershed( inputIP, maskIP, connectivity );
					final long t1 = System.currentTimeMillis();
					ImageStatistics stats = ImageStatistics.getStatistics(resultIP8bit);
					IJ.log(dim + " " + connectivity + " " + "  8 bit " + (t1 - t0) / 1000.0 + " sec mean = " + stats.mean);
					assertEquals( "Different results for 8 bit image of size = " + dim + " (connectivity = " + connectivity + ")",
						0, Double.compare(stats.mean, means.get("[" + dim + "][" + connectivity + "][8]").doubleValue()));
				}
			}

			for( int connectivity : connectivityValues )
			{
				input = copy;

				ImageProcessor inputIP = input.getChannelProcessor();
				
				ImageProcessor maskIP = mask.getChannelProcessor();

				ImageConverter ic = new ImageConverter( input );
				ic.convertToGray16();

				ImageProcessor resultIP16bit;
				{
					final long t0 = System.currentTimeMillis();
					resultIP16bit = Watershed.computeWatershed( inputIP, maskIP, connectivity );
					final long t1 = System.currentTimeMillis();
					ImageStatistics stats = ImageStatistics.getStatistics(resultIP16bit);
					IJ.log(dim + " " + connectivity + " " + " 16 bit " + (t1 - t0) / 1000.0 + " sec mean = " + stats.mean);
					assertEquals( "Different results for 16 bit image of size = " + dim + " (connectivity = " + connectivity + ")",
						0, Double.compare(stats.mean, means.get("[" + dim + "][" + connectivity + "][16]").doubleValue()));
				}
			}

			for( int connectivity : connectivityValues )
			{
				input = copy;

				ImageProcessor inputIP = input.getChannelProcessor();
				
				ImageProcessor maskIP = mask.getChannelProcessor();

				ImageConverter ic = new ImageConverter( input );
				ic.convertToGray32();

				ImageProcessor resultIP32bit;
				{
					final long t0 = System.currentTimeMillis();
					resultIP32bit = Watershed.computeWatershed( inputIP, maskIP, connectivity );
					final long t1 = System.currentTimeMillis();
					ImageStatistics stats = ImageStatistics.getStatistics(resultIP32bit);
					IJ.log(dim + " " + connectivity + " " + " 32 bit " + (t1 - t0) / 1000.0 + " sec mean = " + stats.mean);
					assertEquals( "Different results for32 bit image of size = " + dim + " (connectivity = " + connectivity + ")",
						0, Double.compare(stats.mean, means.get("[" + dim + "][" + connectivity + "][32]").doubleValue()));
				}
			}
		}
	}

	private int diffImagePlus(final ImagePlus a, final ImagePlus b) {
		final int[] dimsA = a.getDimensions(), dimsB = b.getDimensions();
		if (dimsA.length != dimsB.length) return dimsA.length - dimsB.length;
		for (int i = 0; i < dimsA.length; i++) {
			if (dimsA[i] != dimsB[i]) return dimsA[i] - dimsB[i];
		}
		int count = 0;
		final ImageStack stackA = a.getStack(), stackB = b.getStack();
		for (int slice = 1; slice <= stackA.getSize(); slice++) {
			count += diff( stackA.getProcessor( slice ), stackB.getProcessor( slice ) );
		}
		return count;
	}

	private int diff(final ImageProcessor a, final ImageProcessor b) {
		int count = 0;
		final int width = a.getWidth(), height = a.getHeight();
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				if (a.getf(x, y) != b.getf(x, y)) count++;
			}
		}
		return count;
	}
		
	/**
	 * Main method to test and debug the watershed
	 *  
	 * @param args
	 */
	public static void main( final String[] args )
	{
		ImageJ.main( args );
		
		IJ.open( TestWatershed2D.class.getResource( "/files/grains.tif" ).getFile() );
		
		new Watershed3DPlugin().run( null );
	}

}
