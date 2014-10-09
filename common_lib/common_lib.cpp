#include "commonlib_aux.h" 

namespace CommonLib {

//the graphics functions only work in windows
#ifdef WIN32
	#include <windows.h>
	#include <assert.h>

	PBITMAPINFO CreateBitmapInfoStruct(HBITMAP hBmp)
	//auxiliary function for CreateBMPFile
	//don't remember where I copied this from
	{ 
		int errorLevel;
		BITMAP bmp; 
		PBITMAPINFO pbmi; 
		WORD    cClrBits; 

		// Retrieve the bitmap color format, width, and height. 
		errorLevel = GetObject(hBmp, sizeof(BITMAP), (LPSTR)&bmp);
		assert(errorLevel != 0);

		// Convert the color format to a count of bits. 
		cClrBits = (WORD)(bmp.bmPlanes * bmp.bmBitsPixel); 
		if (cClrBits == 1) 
			cClrBits = 1; 
		else if (cClrBits <= 4) 
			cClrBits = 4; 
		else if (cClrBits <= 8) 
			cClrBits = 8; 
		else if (cClrBits <= 16) 
			cClrBits = 16; 
		else if (cClrBits <= 24) 
			cClrBits = 24; 
		else cClrBits = 32; 

		// Allocate memory for the BITMAPINFO structure. (This structure 
		// contains a BITMAPINFOHEADER structure and an array of RGBQUAD 
		// data structures.) 

		 if (cClrBits != 24) 
			 pbmi = (PBITMAPINFO) LocalAlloc(LPTR, 
						sizeof(BITMAPINFOHEADER) + 
						sizeof(RGBQUAD) * (1<< cClrBits)); 

		 // There is no RGBQUAD array for the 24-bit-per-pixel format. 

		 else 
			 pbmi = (PBITMAPINFO) LocalAlloc(LPTR, 
						sizeof(BITMAPINFOHEADER)); 

		// Initialize the fields in the BITMAPINFO structure. 

		pbmi->bmiHeader.biSize = sizeof(BITMAPINFOHEADER); 
		pbmi->bmiHeader.biWidth = bmp.bmWidth; 
		pbmi->bmiHeader.biHeight = bmp.bmHeight; 
		pbmi->bmiHeader.biPlanes = bmp.bmPlanes; 
		pbmi->bmiHeader.biBitCount = bmp.bmBitsPixel; 
		if (cClrBits < 24) 
			pbmi->bmiHeader.biClrUsed = (1<<cClrBits); 

		// If the bitmap is not compressed, set the BI_RGB flag. 
		pbmi->bmiHeader.biCompression = BI_RGB; 

		// Compute the number of bytes in the array of color 
		// indices and store the result in biSizeImage. 
		// For Windows NT, the width must be DWORD aligned unless 
		// the bitmap is RLE compressed. This example shows this. 
		// For Windows 95/98/Me, the width must be WORD aligned unless the 
		// bitmap is RLE compressed.
		pbmi->bmiHeader.biSizeImage = ((pbmi->bmiHeader.biWidth * cClrBits +31) & ~31) /8
									  * pbmi->bmiHeader.biHeight; 
		// Set biClrImportant to 0, indicating that all of the 
		// device colors are important. 
		 pbmi->bmiHeader.biClrImportant = 0; 
		 return pbmi; 
	 } 

	extern "C" void CreateBMPFile(char *pszFile, HBITMAP hBMP, HDC hDC) 
	//given a handle hBMP of a bitmap that is currently selected in the display context hDC,
	//creates a .bmp file named pszFile containing the bitmap
	//don't remember where I copied this from
	{ 
		int errorLevel;
		HANDLE hf;                  // file handle 
		PBITMAPINFO pbi;			// bitmap info
		BITMAPFILEHEADER hdr;       // bitmap file-header 
		PBITMAPINFOHEADER pbih;     // bitmap info-header 
		LPBYTE lpBits;              // memory pointer 
		DWORD dwTotal;              // total count of bytes 
		DWORD cb;                   // incremental count of bytes 
		BYTE *hp;                   // byte pointer 
		DWORD dwTmp; 
		
		pbi = CreateBitmapInfoStruct(hBMP);
		pbih = (PBITMAPINFOHEADER) pbi; 
		lpBits = (LPBYTE) GlobalAlloc(GMEM_FIXED, pbih->biSizeImage);

		assert(lpBits);

		// Retrieve the color table (RGBQUAD array) and the bits 
		// (array of palette indices) from the DIB. 
		errorLevel = GetDIBits(hDC, hBMP, 0, (WORD) pbih->biHeight, lpBits, pbi, DIB_RGB_COLORS);
		assert(errorLevel != 0);

		// Create the .BMP file. 
		hf = CreateFile(pszFile, 
					   GENERIC_READ | GENERIC_WRITE, 
					   (DWORD) 0, 
						NULL, 
					   CREATE_ALWAYS, 
					   FILE_ATTRIBUTE_NORMAL, 
					   (HANDLE) NULL); 
		assert(hf != INVALID_HANDLE_VALUE);
		hdr.bfType = 0x4d42;        // 0x42 = "B" 0x4d = "M" 
		// Compute the size of the entire file. 
		hdr.bfSize = (DWORD) (sizeof(BITMAPFILEHEADER) + 
					 pbih->biSize + pbih->biClrUsed 
					 * sizeof(RGBQUAD) + pbih->biSizeImage); 
		hdr.bfReserved1 = 0; 
		hdr.bfReserved2 = 0; 

		// Compute the offset to the array of color indices. 
		hdr.bfOffBits = (DWORD) sizeof(BITMAPFILEHEADER) + 
						pbih->biSize + pbih->biClrUsed 
						* sizeof (RGBQUAD); 

		// Copy the BITMAPFILEHEADER into the .BMP file. 
		errorLevel = WriteFile(hf, (LPVOID) &hdr, sizeof(BITMAPFILEHEADER), (LPDWORD) &dwTmp,  NULL);
		assert(errorLevel != 0);

		// Copy the BITMAPINFOHEADER and RGBQUAD array into the file. 
		errorLevel = WriteFile(hf, (LPVOID) pbih, sizeof(BITMAPINFOHEADER) 
					  + pbih->biClrUsed * sizeof (RGBQUAD), 
					  (LPDWORD) &dwTmp, NULL);
		assert(errorLevel != 0);

		// Copy the array of color indices into the .BMP file. 
		dwTotal = cb = pbih->biSizeImage; 
		hp = lpBits; 
		errorLevel = WriteFile(hf, (LPSTR) hp, (int) cb, (LPDWORD) &dwTmp,NULL);
		assert(errorLevel != 0);

		// Close the .BMP file. 
		errorLevel = CloseHandle(hf);
		assert(errorLevel != 0);

		// Free memory. 
		GlobalFree((HGLOBAL)lpBits);
	}

	//auxiliary things for HLSToRGB
	const int RGBMAX = 255;
	const int HLSMAX = 240;

	int HueToRGB(int n1, int n2,int hue)
	{

		while (hue > HLSMAX)
			hue -= HLSMAX;
		while (hue < 0)
			hue += HLSMAX;
		if (hue < HLSMAX / 6)
			return ( n1 + (((n2-n1)*hue+(HLSMAX / 12)) / (HLSMAX / 6)));
		else if (hue < HLSMAX / 2)
			return n2;
		else if (hue < (HLSMAX*2)/ 3)
			return n1 + (((n2-n1)*(((HLSMAX*2) / 3)-hue)+(HLSMAX / 12)) / (HLSMAX / 6));
		else 
			return n1;
	}

	COLORREF HLSToRGB (int hue, int lum, int sat)
	//converts a color given as (h,l,s) to RGB(r,g,b)
	{
		int r,g,b;
		int Magic1, Magic2;

		if (sat == 0) {
			r = (lum * RGBMAX) / HLSMAX;
			g = r;
			b = r;
		} else {
			if (lum <= HLSMAX / 2)
				Magic2 = (lum*(HLSMAX + sat) + (HLSMAX / 2)) / HLSMAX;
			else
				Magic2 = lum + sat - ((lum*sat) + (HLSMAX / 2)) / HLSMAX;
			Magic1 = 2*lum-Magic2;
			r = (HueToRGB (Magic1,Magic2,hue+(HLSMAX / 3))*RGBMAX + (HLSMAX / 2)) / HLSMAX;
			g = (HueToRGB (Magic1,Magic2,hue)*RGBMAX + (HLSMAX / 2)) / HLSMAX;
			b = (HueToRGB (Magic1,Magic2,hue-(HLSMAX / 3))*RGBMAX + (HLSMAX / 2)) / HLSMAX;
		}
		return RGB((unsigned char)r, (unsigned char)g, (unsigned char)b);
	}

#endif //WIN32

}; //end namespace CommonLib


