package inra.ijpb.binary.distmap;

import static org.junit.Assert.*;
import ij.ImageStack;
import inra.ijpb.binary.ChamferWeights3D;

import org.junit.Test;

public class DistanceTransform3DShortTest
{

	@Test
	public void testDistanceMap()
	{
		// create 3D image containing a cube 
		ImageStack image = ImageStack.create(20, 20, 20, 8);
		for (int z = 2; z < 19; z++)
		{
			for (int y = 2; y < 19; y++)
			{
				for (int x = 2; x < 19; x++)
				{
					image.setVoxel(x, y, z, 255);
				}
			}
		}
		short[] weights = ChamferWeights3D.BORGEFORS.getShortWeights();
		DistanceTransform3D algo = new DistanceTransform3DShort(weights, true);
		
		ImageStack result = algo.distanceMap(image);
		assertEquals(16, result.getBitDepth());
		
//		System.out.println("result:");
//		for (int x = 0; x < 100; x++)
//		{
//			System.out.print(((int)result.getVoxel(x, 50, 50)) + " ");
//		}
		double middle = result.getVoxel(10, 10, 10);
		assertEquals(9, middle, .1);
	}

	@Test
	public void testDistanceMap_FromCenter()
	{
		// create 3D image containing a cube 
		ImageStack image = ImageStack.create(21, 21, 21, 8);
		for (int z = 0; z < 21; z++)
		{
			for (int y = 0; y < 21; y++)
			{
				for (int x = 0; x < 21; x++)
				{
					image.setVoxel(x, y, z, 255);
				}
			}
		}
		image.setVoxel(10, 10, 10, 0);

		short[] weights = ChamferWeights3D.BORGEFORS.getShortWeights();
		DistanceTransform3D algo = new DistanceTransform3DShort(weights, true);
		
		ImageStack result = algo.distanceMap(image);
		assertEquals(16, result.getBitDepth());
		
		assertEquals(1, result.getVoxel( 9, 10, 10), .1);
		assertEquals(1, result.getVoxel(11, 10, 10), .1);
		assertEquals(4/3, result.getVoxel( 9,  9, 10), .1);
		assertEquals(Math.round(5./3.), result.getVoxel( 9,  9,  9), .1);
		
		// Test some voxels at the cube corners
		int exp = (int) Math.round(10.0 * 5.0 / 3.0);
		assertEquals(exp, result.getVoxel( 0,  0,  0), .01);
		assertEquals(exp, result.getVoxel(20,  0,  0), .01);
		assertEquals(exp, result.getVoxel( 0, 20,  0), .01);
		assertEquals(exp, result.getVoxel(20, 20,  0), .01);
		assertEquals(exp, result.getVoxel( 0,  0, 20), .01);
		assertEquals(exp, result.getVoxel(20,  0, 20), .01);
		assertEquals(exp, result.getVoxel( 0, 20, 20), .01);
		assertEquals(exp, result.getVoxel(20, 20, 20), .01);
	}

}
