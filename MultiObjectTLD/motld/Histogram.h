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

#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <cmath>
#include <cstring>
#include <vector>
#include "Matrix.h"

#define MAX3(a,b,c) a > b ? (a > c ? a : c) : (b > c ? b : c);
#define MIN3(a,b,c) a < b ? (a < c ? a : c) : (b < c ? b : c);

#define NUM_BINS 7
#define GRAY_THRESHOLD 0.125

/// ultra discrete color histograms intended to be used as (weak) classifier asset
class Histogram {
public:
  /// get instance of histogram generating singleton
  static Histogram * getInstance();
  /// creates histogram from whole image
  float * getColorDistribution(const unsigned char * const rgb, const int & size) const;
  /// creates histogram from whole image
  float * getColorDistribution(const unsigned char * const rgb, const int & width, const int & height) const;
  /// creates histogram from image section
  float * getColorDistribution(const unsigned char * const rgb, const int & width, const int & height, const ObjectBox & box) const;
  /// creates a debug image with colors which are maped to same histogram value
  unsigned char * debugImage(const int & bin, int & sideLength) const;
  /// compares two histograms by performing a normalized chi-squared test on their average
  static float compareColorDistribution(const float * const hist1, const float * const hist2);
  
private:
  Histogram();
  ~Histogram();
  
  static void toHS(const float & r, const float & g, const float & b, float & h, float & s);
  static float chiSquareSym(const float * const distr1, const float * const distr2, const int & n);
  static float chiSquare(const float * const correctHistogram, const float * const toCheck, const int & n);
  
  static Histogram * ivInstance;
  unsigned char * ivLookupRGB;
};

Histogram * Histogram::ivInstance = NULL;

Histogram * Histogram::getInstance()
{
  if (ivInstance == NULL)
  {
    ivInstance = new Histogram();
  }
  return ivInstance;
}

Histogram::Histogram()
{
  ivLookupRGB = new unsigned char[4096];
  
  int numColors = 2*(NUM_BINS - 1);
  float colorSize  = 360 / numColors;
  float colorStart = 0.5 * colorSize;
  
#pragma omp parallel for
  for (int c = 0; c < 4096; ++c)
  {
    // determine color voxel
    unsigned char red   = (c >> 8) & 15;
    unsigned char green = (c >> 4) & 15;
    unsigned char blue  =  c       & 15;
    
    // calculate r/g/b value in center of color voxel
    float r = (red   + 0.5) / 16.;
    float g = (green + 0.5) / 16.;
    float b = (blue  + 0.5) / 16.;
    
    // calculate hue and grayness value
    float h, s;
    toHS(r, g, b, h, s);
    
    if (s < GRAY_THRESHOLD)
    {
      ivLookupRGB[c] = 0;
    }
    
    else
    {
      int bin = 1 + (int)((h + colorStart) / colorSize) % numColors;
      ivLookupRGB[c] = bin;
    }
  }  
}

Histogram::~Histogram()
{
  delete [] ivLookupRGB;
  delete ivInstance;
}

float * Histogram::getColorDistribution(const unsigned char * const rgb, const int & size) const
{
  float * result = new float[NUM_BINS]; memset(result, 0., NUM_BINS * sizeof(float));
  
  int numPrebins = 2 * NUM_BINS - 1;
  float * prebinning = new float[numPrebins]; memset(prebinning, 0., numPrebins * sizeof(float));
  
  const unsigned char * red = rgb;
  const unsigned char * green = red + size;
  const unsigned char * blue = green + size;
  
  for (int p = 0; p < size; ++p)
  {
    int r = red[p] >> 4, g = green[p] >> 4, b = blue[p] >> 4;
    int c = (((r<<4)|g)<<4)|b;
    prebinning[ivLookupRGB[c]] += 1;
  }
  
  for (int i = 0; i < numPrebins; ++i)
  {
    prebinning[i] /= size;
  }
  
  // 0  1  2  3  4  5  6  7  8  9 10 11 12
  //    1     2     3     4     5     6
  
  result[0] = prebinning[0];
  prebinning[0] = prebinning[12];
  for (int i = 1; i < NUM_BINS; ++i)
  {
    result[i] = 0.5 * prebinning[2*i-2] + prebinning[2*i-1] + 0.5 * prebinning[2*i];
  }
  
  delete[] prebinning;
  return result;
}

float * Histogram::getColorDistribution(const unsigned char * const rgb, const int & width, const int & height) const
{
  return getColorDistribution(rgb, width * height);
}

float * Histogram::getColorDistribution(const unsigned char * const rgb, const int & width,
                                        const int & height, const ObjectBox & box) const
{
  float * result = new float[NUM_BINS]; memset(result, 0., NUM_BINS * sizeof(float));
  
  int numPrebins = 2 * NUM_BINS - 1;
  float * prebinning = new float[numPrebins]; memset(prebinning, 0., numPrebins * sizeof(float));
  
  int size = width * height;
  const unsigned char* red = rgb;
  const unsigned char* green = red + size;
  const unsigned char* blue = green + size;
  
  int bw = box.width-2, bh = box.height-2;
  int npixels = 0;
  for (int dx = 2; dx < bw; ++dx)
  {
    for (int dy = 2; dy < bh; ++dy)
    {
      int x = round(box.x) + dx, y = round(box.y) + dy;
      if(x < 0 || y < 0)
        continue;
      if(x >= width || y >= height)
        break;
      int p = y * width + x;
      int r = red[p] >> 4, g = green[p] >> 4, b = blue[p] >> 4;
      int c = (((r<<4)|g)<<4)|b;
      prebinning[ivLookupRGB[c]] += 1;
      npixels++;
    }
  }
  
  float invBoxSize = 1.0/npixels;
  for (int i = 0; i < numPrebins; ++i)
  {
    prebinning[i] *= invBoxSize;
  }
  
  result[0] = prebinning[0];
  prebinning[0] = prebinning[12];
  for (int i = 1; i < NUM_BINS; ++i)
  {
    result[i] = 0.5 * prebinning[2*i-2] + prebinning[2*i-1] + 0.5 * prebinning[2*i];
  }
  
  delete[] prebinning;
  return result;
}

float Histogram::compareColorDistribution(const float * const correctHistogram, const float * const candidateHistogram)
{
  return 1 - chiSquareSym(correctHistogram, candidateHistogram, NUM_BINS);
}

unsigned char * Histogram::debugImage(const int & bin, int & sideLength) const
{
  std::vector<unsigned char> rs;
  std::vector<unsigned char> gs;
  std::vector<unsigned char> bs;
  
  for (int c = 0; c < 4096; ++c)
  {
    if (ivLookupRGB[c] == bin)
    {
      // determine color voxel
      unsigned char r = (c >> 8) & 15;
      unsigned char g = (c >> 4) & 15;
      unsigned char b =  c       & 15;
      
      rs.push_back(r * 16 + 0.5);
      gs.push_back(g * 16 + 0.5);
      bs.push_back(b * 16 + 0.5);
    }
  }
  
  int sqsize = ceil(sqrt((float)rs.size()));
  sideLength = sqsize*12;
  
  unsigned char * result = new unsigned char[3*12*12*sqsize*sqsize];
  memset(result,0,3*12*12*sqsize*sqsize*sizeof(unsigned char));
  unsigned char * red   = result;
  unsigned char * green = red   + 12*12*sqsize*sqsize;
  unsigned char * blue  = green + 12*12*sqsize*sqsize;
  
  for (unsigned int i = 0; i < rs.size(); ++i)
  {
    unsigned int row = i / sqsize;
    unsigned int col = i % sqsize;
    
    for (unsigned int y = row * 12 + 1; y < row * 12 + 11; ++y)
    {
      for (unsigned int x = col * 12 + 1; x < col * 12 +11; ++x)
      {
        red  [12 * sqsize * y + x] = rs[i];
        green[12 * sqsize * y + x] = gs[i];
        blue [12 * sqsize * y + x] = bs[i];
      }
    }
  }
  return result;
}

void Histogram::toHS(const float & r, const float & g, const float & b, float & h, float & s)
{
  float max = MAX3(r,g,b);
  float min = MIN3(r,g,b);
  float dmm = max - min;
  
  if      (dmm == 0) h = 0; 
  else if (max == r) h = 60 * (0 + (g-b)/dmm);
  else if (max == g) h = 60 * (2 + (b-r)/dmm);
  else if (max == b) h = 60 * (4 + (r-g)/dmm);
  h = h < 0 ? h + 360 : h;
  
  if (max == 0) s = 0;
  else          s = dmm / max;
}

inline float Histogram::chiSquareSym(const float * const v1, const float * const v2, const int & n)
{
  float chisqr2 = 0.;
  for (int j = 0; j < n; ++j)
  {
    float diff = v1[j] - v2[j];
    float sum  = v1[j] + v2[j];
    chisqr2 += sum == 0 ? 0 : diff * diff / sum;
  }
  return 0.5 * chisqr2;
}

float Histogram::chiSquare(const float*const correctHistogram, const float*const toCheck, const int& n)
{
  float chisqr = 0.;
  for (int j = 0; j < n; ++j)
  {
    float correct = correctHistogram[j];
    float diff = toCheck[j] - correct;
    chisqr += correct > 0.0001 ? diff * diff / correct : 0;
  }
  return chisqr;
  
}


#endif // HISTOGRAM_H
