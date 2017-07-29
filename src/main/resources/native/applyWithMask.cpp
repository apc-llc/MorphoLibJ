#include <algorithm>
#include <cstdio>
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
		if (value != other.value)
			return value < other.value;
		
		if (x != other.x)
			return x < other.x;
		
		if (y != other.y)
			return y < other.y;
	}
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
	jobject pixelListObj, jobjectArray tabLabelsObj, jobject maskImageObj,
	jint MASK, jint WSHED, jint INIT, jint INQUEUE)
{
    int currentLabel = 0;
    
    bool flag = false;	    

	vector<Pixel> pixels;
	{
		if( verbose ) IJ.log("  Converting pixels list to native array..." );
		long t1 = System.currentTimeMillis();
		IJ.showStatus("Converting pixels list to native array...");

		jclass arrayClass = env->FindClass("java/util/ArrayList");
		jmethodID sizeMethodID = env->GetMethodID(arrayClass, "size", "()I");
		jclass integerClass = env->FindClass("java/lang/Integer");
		jmethodID intValueMethodID = env->GetMethodID(integerClass, "intValue", "()I");
		jint size = env->CallIntMethod(pixelListObj, sizeMethodID);

		pixels.resize(size);

		jmethodID getMethodID = env->GetMethodID(arrayClass, "get", "(I)Ljava/lang/Object;");
		jclass pixelRecordClass = env->FindClass("inra/ijpb/data/PixelRecord");
		jfieldID valueFieldID = env->GetFieldID(pixelRecordClass, "value", "D");
		jfieldID cursor2DFieldID = env->GetFieldID(pixelRecordClass, "cursor", "Linra/ijpb/data/Cursor2D;");
		jclass cursor2DClass = env->FindClass("inra/ijpb/data/Cursor2D");
		jmethodID getXMethodID = env->GetMethodID(cursor2DClass, "getX", "()I");
		jmethodID getYMethodID = env->GetMethodID(cursor2DClass, "getY", "()I");
		for (int i = 0; i < size; i++)
		{
			jobject element = env->CallObjectMethod(pixelListObj, getMethodID, i);
			
			Pixel& pixel = pixels[i];
			pixel.value = env->GetDoubleField(element, valueFieldID);

			jobject cursor;
			cursor = env->GetObjectField(element, cursor2DFieldID);
			pixel.x = env->CallIntMethod(cursor, getXMethodID);
			pixel.y = env->CallIntMethod(cursor, getYMethodID);			

			env->DeleteLocalRef(element);
		}

	    long t2 = System.currentTimeMillis();
	    stringstream ss;
	    ss << "  Converting took " << (t2-t1) << " ms.";
	    if( verbose ) IJ.log(ss.str());
	}
	
	vector<float> maskImage(size1 * size2);
	{
		if( verbose ) IJ.log("  Converting mask image to native array..." );
		long t1 = System.currentTimeMillis();
		IJ.showStatus("Converting mask image to native array...");

		jclass imageProcessorClass = env->FindClass("ij/process/ImageProcessor");
		jmethodID getfMethodID = env->GetMethodID(imageProcessorClass, "getf", "(II)F"); 

		for (int i = 0; i < size1; i++)
			for (int j = 0; j < size2; j++)
				maskImage[i * size2 + j] = env->CallFloatMethod(maskImageObj, getfMethodID, i, j);

	    long t2 = System.currentTimeMillis();
	    stringstream ss;
	    ss << "  Converting took " << (t2-t1) << " ms.";
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
	vector<int> tabLabels(size1 * size2, INIT);
	{
		IJ.log( "  Flooding..." );
		IJ.showStatus( "Flooding..." );
		long start = System.currentTimeMillis();

		int currentIndex = 0;
		int heightIndex1 = currentIndex;
		int heightIndex2 = currentIndex;
		
		queue<Cursor2D> fifo;

		vector<Cursor2D> neighs(8);
		if (connectivity >= 4)
		{
			neighs[0] = Cursor2D(-1, 0);
			neighs[1] = Cursor2D(1, 0);
			neighs[2] = Cursor2D(0, -1);
			neighs[3] = Cursor2D(0, 1); 
		}
		if (connectivity == 8)
		{
			neighs[4] = Cursor2D(1, 1);
			neighs[5] = Cursor2D(-1, 1);
			neighs[6] = Cursor2D(1, -1);
			neighs[7] = Cursor2D(-1, -1); 
		}

		// for h <- h_min to h_max; geodesic SKIZ of level h-1 inside level h
		while( currentIndex < pixels.size() )
		{
			double h = pixels[currentIndex].value;

			for(int pixelIndex = heightIndex1, e = pixels.size(); pixelIndex < e; pixelIndex ++)
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
				tabLabels[ i * size2 + j ] = MASK;

				// read neighbor coordinates
				for (int k = 0; k < connectivity; k++)
				{
					int u = i + neighs[k].x;
					int v = j + neighs[k].y;

					// initialize queue with neighbors at level h of current basins or watersheds
					if ( u >= 0 && u < size1 && v >= 0 && v < size2
						&& tabLabels[ u * size2 + v ] >= WSHED && maskImage[ u * size2 + v ] > 0)
					//	&&  ( tabLabels[ u ][ v ] > 0 || tabLabels[ u ][ v ] == WSHED ) )
					{
						fifo.push(Cursor2D(i, j));
						tabLabels[ i * size2 + j ] = INQUEUE;
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

					if ( u >= 0 && u < size1 && v >= 0 && v < size2 && maskImage[ u * size2 + v ] > 0)
					{
						if ( tabLabels[ u * size2 + v ] > 0 ) // i.e. the pixel belongs to an already labeled basin
						{
							if ( tabLabels[ i * size2 + j ] == INQUEUE || (tabLabels[ i * size2 + j ] == WSHED && flag == true ) )
							{
								tabLabels[ i * size2 + j ] = tabLabels[ u * size2 + v ];
							}
							else if ( tabLabels[ i * size2 + j ] > 0 && tabLabels[ i * size2 + j ] != tabLabels[ u * size2 + v ] )
							{
								tabLabels[ i * size2 + j ] = WSHED;
								flag = false;
							}
						}
						else if ( tabLabels[ u * size2 + v ] == WSHED )
						{
							if( tabLabels[ i * size2 + j ] == INQUEUE )
							{
								tabLabels[ i * size2 + j ] = WSHED;
								flag = true;
							}
						}
						else if ( tabLabels[ u * size2 + v ] == MASK )
						{
							tabLabels[ u * size2 + v ] = INQUEUE;
							fifo.push(Cursor2D(u, v));
						}
					}
				}
			}

			// check for new minima at level h
			
			for(int pixelIndex = heightIndex2, e = pixels.size(); pixelIndex < e; pixelIndex ++, currentIndex++)
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

				if ( tabLabels[ i * size2 + j ] == MASK ) // the pixel is inside a new minimum
				{
					currentLabel ++;
					fifo.push(Cursor2D(i, j));
					tabLabels[ i * size2 + j ] = currentLabel;
					
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
								&& tabLabels[ u * size2 + v ] == MASK && maskImage[ u * size2 + v ] > 0)
							{
								fifo.push(Cursor2D(u, v));
								tabLabels[ u * size2 + v ] = currentLabel;
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
		if( verbose ) IJ.log("  Converting tab labels to managed array..." );
		long t1 = System.currentTimeMillis();
		IJ.showStatus("Converting tab labels to managed array...");

		for (int i = 0; i < size1; i++)
		{
			jintArray jline = env->NewIntArray(size2);
			env->SetIntArrayRegion(jline, 0, size2, &tabLabels[i * size2]);
			env->SetObjectArrayElement(tabLabelsObj, i, jline);
			env->DeleteLocalRef(jline);
		}

	    long t2 = System.currentTimeMillis();
	    stringstream ss;
	    ss << "  Converting took " << (t2-t1) << " ms.";
	    if( verbose ) IJ.log(ss.str());
	}
}

