/* Copyright (C) 2012 Christian Lutz, Thorsten Engesser
 * 
 * This file is part of motld
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef UTILS_H
#define UTILS_H

#include <cstring>

#ifdef _MSC_VER
#include <time.h>
#else
#include <sys/time.h>
#endif

/*****************************************************************************
 *                                 General Stuff                             *
 *****************************************************************************/

/// returns random integer n with min <= n <= max
inline int randInt(int min, int max)
{
  return min + rand() % (1 + max - min);
}

/// returns random float x with min <= x <= max
inline float randFloat(float min, float max)
{
  return min + ((float)rand() / RAND_MAX)*(max-min);
}

/// returns time value in milliseconds
inline int getTime()
{
  #ifdef _MSC_VER
  clock_t t = clock();
  return t * 1000 / CLOCKS_PER_SEC;
  #else
  timeval t;
  gettimeofday(&t, 0);
  return t.tv_sec*1000 + (int)(t.tv_usec/1000.);
  #endif
}

/*****************************************************************************
 *                   File Input / Output - PGM and PPM                       *
 *****************************************************************************/

/// reads from PPM file into array of something
template <class T>
T* readFromPPM(const char* aFilename, int & xSize, int & ySize, int & zSize) {
  T * result = 0;
  FILE *aStream;
  aStream = fopen(aFilename,"rb");
  if (aStream == 0)
    std::cerr << "File not found: " << aFilename << std::endl;
  int dummy;
  // Find beginning of file (P6)
  while (getc(aStream) != 'P');
  dummy = getc(aStream);
  if (dummy == '5') zSize = 1;
  else if (dummy == '6') zSize = 3;
  else
  {
    std::cerr << "Cannot read File - Invalid File Format!" << std::endl;
    zSize = 0;
    return result;
  }
  do dummy = getc(aStream); while (dummy != '\n' && dummy != ' ');
  // Remove comments and empty lines
  dummy = getc(aStream);
  while (dummy == '#') {
    while (getc(aStream) != '\n');
    dummy = getc(aStream);
  }
  while (dummy == '\n')
    dummy = getc(aStream);
  // Read image size
  xSize = dummy-48;
  while ((dummy = getc(aStream)) >= 48 && dummy < 58)
    xSize = 10*xSize+dummy-48;
  while ((dummy = getc(aStream)) < 48 || dummy >= 58);
  ySize = dummy-48;
  while ((dummy = getc(aStream)) >= 48 && dummy < 58)
    ySize = 10*ySize+dummy-48;
  while (dummy != '\n' && dummy != ' ')
    dummy = getc(aStream);
  while (dummy < 48 || dummy >= 58) dummy = getc(aStream);
  while ((dummy = getc(aStream)) >= 48 && dummy < 58);
  if (dummy != '\n') while (getc(aStream) != '\n');
  // Adjust size of data structure
  result = new T[xSize*ySize*zSize];
  //result = (T*)malloc(xSize*ySize*zSize*sizeof(T));
  // Read image data
  int aSize = xSize*ySize;
  if (zSize == 1)
    for (int i = 0; i < aSize; i++)
      result[i] = getc(aStream) / 255.;
  else {
    int aSizefloatwice = aSize+aSize;
    for (int i = 0; i < aSize; i++) {
      result[i] = getc(aStream);
      result[i+aSize] = getc(aStream);
      result[i+aSizefloatwice] = getc(aStream);
    }
  }
  fclose(aStream);
  return result;
}

/// reads from PGM file into array of something
template<class T>
T* readFromPGM(const char* filename, int & xSize, int & ySize)
{
  T * result = 0;
  FILE *flStream;
  flStream = fopen(filename,"rb");
  if (flStream == 0)
  { 
    std::cerr << "File not found: " << filename << std::endl;
    return result;
  }
  int dummy;
  // Find beginning of file (P5)
  while (getc(flStream) != 'P');
  if (getc(flStream) != '5')
  {
    std::cerr << "Cannot read File - Invalid File Format!" << std::endl;
    return result;
  }
  do dummy = getc(flStream); while (dummy != '\n' && dummy != ' ');
  // Remove comments and empty lines
  dummy = getc(flStream);
  while (dummy == '#') {
    while (getc(flStream) != '\n');
    dummy = getc(flStream);
  }
  while (dummy == '\n')
    dummy = getc(flStream);
  // Read image size
  xSize = dummy-48;
  while ((dummy = getc(flStream)) >= 48 && dummy < 58)
    xSize = 10*xSize+dummy-48;
  while ((dummy = getc(flStream)) < 48 || dummy >= 58);
  ySize = dummy-48;
  while ((dummy = getc(flStream)) >= 48 && dummy < 58)
    ySize = 10*ySize+dummy-48;
  while (dummy != '\n' && dummy != ' ')
    dummy = getc(flStream);
  while ((dummy = getc(flStream)) >= 48 && dummy < 58);
  if (dummy != '\n') while (getc(flStream) != '\n');
  result = new T[xSize*ySize];
  for (int i = 0; i < xSize*ySize; i++)
    result[i] = getc(flStream);
  fclose(flStream);
  return result;
}

/// writes from array of something into PGM file
template <class T>
void writeToPGM(const char *filename, T *data, int xSize, int ySize) {
  FILE *flStream;
  flStream = fopen(filename,"wb");
  // write header
  char line[60];
  sprintf(line,"P5\n%d %d\n255\n",xSize,ySize);
  fwrite(line,strlen(line),1,flStream);
  // write data
  for (int i = 0; i < xSize*ySize; i++) {
    char dummy = (char)data[i];
    fwrite(&dummy,1,1,flStream);
  }
  fclose(flStream);
}

/// writes from array of something into PPM file
template <class T>
void writeToPPM(const char * aFilename, T * data, int xSize, int ySize) {
  
  T * red = data;
  T * green = red + xSize * ySize;
  T * blue = green + xSize * ySize;
  
  FILE* outimage = fopen(aFilename, "wb");
  fprintf(outimage, "P6 \n");
  fprintf(outimage, "%d %d \n255\n", xSize,ySize);
  for (int p = 0; p < xSize * ySize; ++p)
  {
    fwrite (red+p, sizeof(unsigned char), 1, outimage);
    fwrite (green+p, sizeof(unsigned char), 1, outimage);
    fwrite (blue+p, sizeof(unsigned char), 1, outimage);
   }
    
  fclose(outimage);
}

/// reduces array with 3 sequential color channels to 1 channel array
template <class T>
T* toGray(T * rgb, int size)
{
  T* result = new T[size];
  T* red = rgb;
  T* green = red + size;
  T* blue = green + size;
  for (int i = 0; i < size; ++i)
  {
    result[i] = (red[i] + green[i] + blue[i]) / 3.;
  }
  return result;
}

#endif // UTILS_H
