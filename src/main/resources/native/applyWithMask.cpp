#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <jni.h>
#include <queue>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

static class
{
public :
	static void showStatus(string msg) { /* TODO */ }
	
	static void log(string msg) { /* TODO */ }
	
	static void showProgress(double progress) { /* TODO */ }
}
IJ;

static class
{
public :
	static long currentTimeMillis() { return 0; }	
}
System;

struct Pixel
{
	int x, y;
	double value;
	
	bool operator<(const Pixel& other)
	{
		// Special case: https://docs.oracle.com/javase/7/docs/api/java/lang/Double.html#equals(java.lang.Object)
		if (signbit(value) && !signbit(other.value))
			return true;

		if (value != other.value)
			return value < other.value;
		
		if (x != other.x)
			return x < other.x;
		
		if (y != other.y)
			return y < other.y;
	}
	
	Pixel() : x(0), y(0), value(0.0) { }
	
	Pixel(int x_, int y_, double value_) : x(x_), y(y_), value(value_) { }
};

struct Cursor2D
{
	int x, y;

	Cursor2D() : x(0), y(0) { }
	
	Cursor2D(int x_, int y_) : x(x_), y(y_) { }
};

extern "C" JNIEXPORT void JNICALL Java_inra_ijpb_watershed_WatershedTransform2D_applyWithMask(
	JNIEnv *env, jclass object, jdouble hMin, jdouble hMax,
	jint size1, jint size2, jint connectivity, jboolean verbose,
	jobjectArray imagePixelsObj, jobjectArray maskPixelsObj, jfloatArray resultPixelsObj,
	jint MASK, jint WSHED, jint INIT, jint INQUEUE)
{
    int currentLabel = 0;
    
    bool flag = false;

	vector<jfloatArray> rowObjs(size2);
	vector<float*> maskPixels(size2);
	{
		if( verbose ) IJ.log("  Converting image mask to native array..." );
		long t1 = System.currentTimeMillis();
		IJ.showStatus("Converting image mask to native array...");

		for (int i = 0; i < size2; i++)
		{
			jfloatArray rowObj = (jfloatArray)env->GetObjectArrayElement(maskPixelsObj, i);
			maskPixels[i] = (float*)env->GetFloatArrayElements(rowObj, 0);
			rowObjs[i] = rowObj;
		}

	    long t2 = System.currentTimeMillis();
	    stringstream ss;
	    ss << "  Converting took " << (t2-t1) << " ms.";
	    if( verbose ) IJ.log(ss.str());
	}

	vector<Pixel> pixels(size1 * size2);
	{
		vector<float> imagePixels(size1 * size2);
		{
			if( verbose ) IJ.log("  Converting image pixels to native array..." );
			long t1 = System.currentTimeMillis();
			IJ.showStatus("Converting image pixels to native array...");

			for (int i = 0; i < size2; i++)
			{
				jfloatArray rowObj = (jfloatArray)env->GetObjectArrayElement(imagePixelsObj, i);
				jfloat* row = env->GetFloatArrayElements(rowObj, 0);
				memcpy(&imagePixels[i * size1], row, sizeof(float) * size1);

				env->DeleteLocalRef(rowObj);
			}

			long t2 = System.currentTimeMillis();
			stringstream ss;
			ss << "  Converting took " << (t2-t1) << " ms.";
			if( verbose ) IJ.log(ss.str());
		}

		if( verbose ) IJ.log("  Filtering masked pixels..." );
		IJ.showStatus("Sorting pixels by value...");
		long t1 = System.currentTimeMillis();
	
		int nfiltered = 0;
		
		for (int j = 0; j < size2; j++)
			for (int i = 0; i < size1; i++)
			{
				double h = imagePixels[j * size1 + i];
				if ((maskPixels[j][i] > 0) && (h >= hMin) && (h <= hMax))
					pixels[nfiltered++] = Pixel(j, i, imagePixels[j * size1 + i]);
			}
		
		pixels.resize(nfiltered);

		long t2 = System.currentTimeMillis();
		stringstream ss;
		ss << "  Filtering took " << (t2-t1) << " ms.";
		if( verbose ) IJ.log(ss.str());
	}

	{
		if( verbose ) IJ.log("  Sorting pixels by value..." );
		IJ.showStatus("Sorting pixels by value...");
		long t1 = System.currentTimeMillis();
		sort(pixels.begin(), pixels.end());
		long t2 = System.currentTimeMillis();
		stringstream ss;
		ss << "  Sorting took " << (t2-t1) << " ms.";
		if( verbose ) IJ.log(ss.str());
    }
    
	// value INIT is assigned to each pixel of the output labels
	int* tabLabels = (int*)env->GetFloatArrayElements(resultPixelsObj, 0);
	{
		for (int i = 0, e = size1 * size2; i < e; i++)
			tabLabels[i] = INIT;
	}

	{
		IJ.log( "  Flooding..." );
		IJ.showStatus( "Flooding..." );
		long start = System.currentTimeMillis();

		int currentIndex = 0;
		int heightIndex1 = currentIndex;
		int heightIndex2 = currentIndex;
		
		queue<Cursor2D> fifo;

		vector<Cursor2D> neighs(8);
		if (connectivity == 4)
		{
			neighs[0] = Cursor2D(-1,  0);
			neighs[1] = Cursor2D( 0, -1);
			neighs[2] = Cursor2D( 1,  0);
			neighs[3] = Cursor2D( 0,  1);
		}
		else if (connectivity == 8)
		{
			neighs[0] = Cursor2D(-1, -1);
			neighs[1] = Cursor2D(-1,  0);
			neighs[2] = Cursor2D(-1,  1);
			neighs[3] = Cursor2D( 0, -1);
			neighs[4] = Cursor2D( 0,  1);
			neighs[5] = Cursor2D( 1, -1);
			neighs[6] = Cursor2D( 1,  0);
			neighs[7] = Cursor2D( 1,  1);
		}

		// for h <- h_min to h_max; geodesic SKIZ of level h-1 inside level h
		while( currentIndex < pixels.size() )
		{
			double h = pixels[currentIndex].value;

			for(int pixelIndex = heightIndex1, e = pixels.size(); pixelIndex < e; pixelIndex++)
			{
				const Pixel& pixel = pixels[pixelIndex];

				if( pixel.value != h )
				{
					// this pixel is at level h+1
					heightIndex1 = pixelIndex;
					break;
				}

				const int& i = pixel.x;
				const int& j = pixel.y;
				
				// set label to MASK
				tabLabels[ j * size1 + i ] = MASK;

				// read neighbor coordinates
				for (int k = 0; k < connectivity; k++)
				{
					int u = i + neighs[k].x;
					int v = j + neighs[k].y;

					// initialize queue with neighbors at level h of current basins or watersheds
					if ( u >= 0 && u < size1 && v >= 0 && v < size2
						&& tabLabels[ v * size1 + u ] >= WSHED && maskPixels[u][v] > 0)
					//	&&  ( tabLabels[ v * size1 + u ] > 0 || tabLabels[ v * size1 + u ] == WSHED ) )
					{
						fifo.push(Cursor2D(i, j));
						tabLabels[ j * size1 + i ] = INQUEUE;
						break;
					}
				}// end for
			}// end for

			while( fifo.empty() == false )
			{
				// retrieve point p
				const Cursor2D p = fifo.front();
				fifo.pop();
				const int& i = p.x;
				const int& j = p.y;

				// read neighbor coordinates
				for(int k = 0; k < connectivity; k++)
				{
					// labeling current point by inspecting neighbors
					int u = i + neighs[k].x;
					int v = j + neighs[k].y;

					if ( u >= 0 && u < size1 && v >= 0 && v < size2 && maskPixels[u][v] > 0)
					{
						if ( tabLabels[ v * size1 + u ] > 0 ) // i.e. the pixel belongs to an already labeled basin
						{
							if ( tabLabels[ j * size1 + i ] == INQUEUE || (tabLabels[ j * size1 + i ] == WSHED && flag == true ) )
							{
								tabLabels[ j * size1 + i ] = tabLabels[ v * size1 + u ];
							}
							else if ( tabLabels[ j * size1 + i ] > 0 && tabLabels[ j * size1 + i ] != tabLabels[ v * size1 + u ] )
							{
								tabLabels[ j * size1 + i ] = WSHED;
								flag = false;
							}
						}
						else if ( tabLabels[ v * size1 + u ] == WSHED )
						{
							if( tabLabels[ j * size1 + i ] == INQUEUE )
							{
								tabLabels[ j * size1 + i ] = WSHED;
								flag = true;
							}
						}
						else if ( tabLabels[ v * size1 + u ] == MASK )
						{
							tabLabels[ v * size1 + u ] = INQUEUE;
							fifo.push(Cursor2D(u, v));
						}
					}
				}
			}

			// check for new minima at level h
			
			for(int pixelIndex = heightIndex2, e = pixels.size(); pixelIndex < e; pixelIndex++, currentIndex++)
			{
				const Pixel& pixel = pixels[pixelIndex];
				
				if( pixel.value != h )
				{
					// this pixel is at level h+1
					heightIndex2 = pixelIndex;
					break;
				}
    			
				const int& i = pixel.x;
				const int& j = pixel.y;

				if ( tabLabels[ j * size1 + i ] == MASK ) // the pixel is inside a new minimum
				{
					currentLabel++;
					fifo.push(Cursor2D(i, j));
					tabLabels[ j * size1 + i ] = currentLabel;
					
					while( fifo.empty() == false )
					{
						const Cursor2D p2 = fifo.front();
						fifo.pop();

						// read neighbor coordinates
						for (int k = 0; k < connectivity; k++)
						{
							int u = p2.x + neighs[k].x;
							int v = p2.y + neighs[k].y;

							if ( u >= 0 && u < size1 && v >= 0 && v < size2
								&& tabLabels[ v * size1 + u ] == MASK && maskPixels[u][v] > 0)
							{
								fifo.push(Cursor2D(u, v));
								tabLabels[ v * size1 + u ] = currentLabel;
							}
						}// end for
					}// end while
				}// end if
			}// end for

			IJ.showProgress( h / hMax );

		}// end while (flooding)

		IJ.showProgress( 1.0 );

		long end = System.currentTimeMillis();
		stringstream ss;
		ss << "  Flooding took: " << (end-start) << " ms";
		if( verbose ) IJ.log(ss.str());
	}

	{
		for (int i = 0; i < size2; i++)
			env->DeleteLocalRef(rowObjs[i]);
	}

	{
		if( verbose ) IJ.log("  Converting tab labels to resulting image pixels..." );
		long t1 = System.currentTimeMillis();
		IJ.showStatus("Converting tab labels to resulting image pixels...");

		// For speed, the resulting pixels array is mapped to the same storage as used for tabLabels.
		// Of course, this assumes sizeof(int) = sizeof(float).
		float* resultPixels = (float*)tabLabels;

		for (int j = 0; j < size2; j++)
		{
			for (int i = 0; i < size1; i++)
			{
				float* resultPixel = &resultPixels[ j * size1 + i ];
				int label = *(int*)resultPixel;

				if ( label == INIT) // set unlabeled pixels to 0
					*resultPixel = 0;
				else
					*resultPixel = label;
			}
		}
		
		env->ReleaseFloatArrayElements(resultPixelsObj, resultPixels, 0);

	    long t2 = System.currentTimeMillis();
	    stringstream ss;
	    ss << "  Converting took " << (t2-t1) << " ms.";
	    if( verbose ) IJ.log(ss.str());
	}
}

