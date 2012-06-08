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

#ifndef NNCLASSIFIER_H
#define NNCLASSIFIER_H

//#include <iostream>
//#include <fstream>
#include <vector>
#include "Matrix.h"
#include "Histogram.h"

/// Data structure representing nearest neighbor patches with color histograms
class NNPatch
{
public:
  /// The patch (will be normalized to zero mean)
  Matrix patch;
  /// The average of the grayscale values (to enable reconstruction)
  float avg;
  /// Saves the squared norm of the grayscale values to speed up computation later
  float norm2;
  /// A normalized color histogram or NULL if not used
  float* histogram;
  /// Default (empty) constructor
  NNPatch() : patch(), avg(0), norm2(1), histogram(NULL){}
  /// Copy constructor
  NNPatch(const NNPatch& copyFrom);
  /// Constructor providing only a patch
  NNPatch(const Matrix& curPatch);  
  /// Constructor creating the color histogram
  NNPatch(const Matrix& curPatch, const ObjectBox& bbox, 
          const unsigned char * rgb = NULL, const int w = 0, const int h = 0);
  /// Constructor extracting the patch and the color histogram out of the image
  NNPatch(const ObjectBox& bbox, const Matrix& curImage, const int patchSize,
          const unsigned char * rgb = NULL, const int w = 0, const int h = 0);
  /// Constructor for loading from file
  NNPatch(std::ifstream & inputStream, const int patchSize);
  /// Destructor
  ~NNPatch();
  /// Copy operator
  NNPatch& operator=(const NNPatch& copyFrom);
  /// Method for saving to file
  void saveToStream(std::ofstream & outputStream) const;
};

/** @brief The nearest neighbor classifier is invoked at the top level to evaluate detections.
 */
class NNClassifier
{
public:
  /// Constructor
  NNClassifier(int width, int height, int patchSize, bool useColor = true, bool allowFastChange = false);
  /// Constructor for loading from file
  NNClassifier(std::ifstream & inputStream);
  /// Returns the confidence of a given patch with respect to a certain class.
  double getConf(const NNPatch& patch, int objId = 0, bool conservative = false) const;
  /// Returns the confidence of a given patch while subsequently computing and saving the color histogram if needed.
  double getConf(NNPatch& patch, int objId, bool conservative,
                  const ObjectBox& bbox, const unsigned char * rgb, int w, int h) const;
  /// Trains a new patch to the classifier if it is considered "new" enough.
  bool trainNN(const NNPatch& patch, int objId = 0, bool positive = true, bool tmp = false);
  /// Initializes a new object class with the given patch.
  void addObject(const NNPatch& patch);
  /// Returns a pointer to positive patches (intended for drawing).
  const std::vector<std::vector<NNPatch> > * getPosPatches() const;
  /// Returns a pointer to negative patches (intended for drawing).
  const std::vector<NNPatch> * getNegPatches() const;
  /// Removes previously added warps (rotated patches) from positive list.
  void removeWarps();
  /// Saves the classifier (i.e. the patches) to file.
  void saveToStream(std::ofstream & outputStream) const;
  
private:
  int ivWidth;
  int ivHeight;
  int ivPatchSize;
  static Histogram * ivHistogram;
  std::vector<std::vector<NNPatch> > ivPosPatches;
  std::vector<NNPatch> ivNegPatches;
  std::vector<char> ivWarpIndices;
  bool ivUseColor, ivAllowFastChange;
  double getConf(const float* patch, float norm2 = 1.0f, int objId = 0, bool conservative = false) const;
  double crossCorr(const float* patchA, const float* patchB, float denom = 1) const;
  double cmpHistograms(const float* h1, const float* h2) const;
};


/**************************************************************************************************
 * IMPLEMENTATION                                                                                 *
 **************************************************************************************************/
 
////////////////////////////////////////////////////////////////////////////////////////////////////
// NNPatch
 
NNPatch::NNPatch(const NNPatch& copyFrom)
{
  patch = Matrix(copyFrom.patch);
  avg = copyFrom.avg;
  norm2 = copyFrom.norm2;
  if(copyFrom.histogram != NULL)
  {
    histogram = new float[NUM_BINS];
    memcpy(histogram, copyFrom.histogram, NUM_BINS * sizeof(float));
    //for(int i = 0; i < NUM_BINS; i++)
    //  histogram[i] = copyFrom.histogram[i];
  }else
    histogram = NULL;
}
  
NNPatch::NNPatch(const Matrix& curPatch)
{
  patch = curPatch;
  avg = patch.avg();
  patch += -avg;
  norm2 = patch.norm2();
  histogram = NULL;
}  

NNPatch::NNPatch(const Matrix& curPatch, const ObjectBox& bbox, 
                 const unsigned char * rgb, const int w, const int h)
{
  patch = curPatch;
  avg = patch.avg();
  patch += -avg;
  norm2 = patch.norm2();
  histogram = (rgb == NULL ? NULL 
                  : Histogram::getInstance()->getColorDistribution(rgb, w, h, bbox));
}

NNPatch::NNPatch(const ObjectBox& bbox, const Matrix& curImage, const int patchSize, 
                 const unsigned char * rgb, const int w, const int h)
{
  patch = curImage.getRectSubPix(bbox.x + 0.5 * bbox.width, bbox.y + 0.5 * bbox.height, 
                                  round(bbox.width), round(bbox.height));
  patch.rescale(patchSize, patchSize);
  avg = patch.avg();
  patch += -avg;
  norm2 = patch.norm2();
  histogram = (rgb == NULL ? NULL 
                  : Histogram::getInstance()->getColorDistribution(rgb, w, h, bbox));
}

NNPatch::NNPatch(std::ifstream & inputStream, const int patchSize)
{
  patch = Matrix(patchSize, patchSize);
  inputStream.read((char*)patch.data(), patchSize*patchSize*sizeof(float));  
  inputStream.read((char*)&avg, sizeof(float));  
  inputStream.read((char*)&norm2, sizeof(float));  
  bool hist;
  inputStream.read((char*)&hist, sizeof(bool));  
  if(hist){
    histogram = new float[NUM_BINS];
    inputStream.read((char*)histogram, NUM_BINS*sizeof(float));  
  }else
    histogram = NULL;
}

void NNPatch::saveToStream(std::ofstream & outputStream) const
{
  outputStream.write((char*)patch.data(), patch.size()*sizeof(float));  
  outputStream.write((char*)&avg, sizeof(float));  
  outputStream.write((char*)&norm2, sizeof(float));  
  bool hist = histogram != NULL;
  outputStream.write((char*)&hist, sizeof(bool));  
  if(hist)
    outputStream.write((char*)histogram, NUM_BINS*sizeof(float));  
}

NNPatch::~NNPatch()
{
  if(histogram != NULL)
    delete [] histogram;
  histogram = NULL;
}

NNPatch& NNPatch::operator=(const NNPatch& copyFrom)
{
  if (this != &copyFrom) {
    patch = Matrix(copyFrom.patch);
    avg = copyFrom.avg;
    norm2 = copyFrom.norm2;
    if(copyFrom.histogram != NULL)
    {
      if(histogram != NULL)
        delete [] histogram;
      histogram = new float[NUM_BINS];
      memcpy(histogram, copyFrom.histogram, NUM_BINS * sizeof(float));
    }else
      histogram = NULL;
  }
  return *this;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// NNClassifier

Histogram * NNClassifier::ivHistogram = Histogram::getInstance();

NNClassifier::NNClassifier(int width, int height, int patchSize, bool useColor, bool allowFastChange) 
  : ivWidth(width), ivHeight(height), ivPatchSize(patchSize),
    ivUseColor(useColor), ivAllowFastChange(allowFastChange) {}

NNClassifier::NNClassifier(std::ifstream & inputStream)
{
  inputStream.read((char*)&ivWidth, sizeof(int));
  inputStream.read((char*)&ivHeight, sizeof(int));
  inputStream.read((char*)&ivPatchSize, sizeof(int));
  inputStream.read((char*)&ivUseColor, sizeof(bool));
  inputStream.read((char*)&ivAllowFastChange, sizeof(bool));
  int nNeg, nObs, nPos;
  inputStream.read((char*)&nNeg, sizeof(int));
  for(int i = 0; i < nNeg; ++i)
    ivNegPatches.push_back(NNPatch(inputStream, ivPatchSize));
  inputStream.read((char*)&nObs, sizeof(int));
  ivPosPatches = std::vector<std::vector<NNPatch> >(nObs);
  for(int i = 0; i < nObs; ++i)
  {
    inputStream.read((char*)&nPos, sizeof(int));
    for(int j = 0; j < nPos; ++j)
      ivPosPatches[i].push_back(NNPatch(inputStream, ivPatchSize));  
  }
}

const std::vector<std::vector<NNPatch> > * NNClassifier::getPosPatches() const
{
  return &ivPosPatches;
}
const std::vector<NNPatch> * NNClassifier::getNegPatches() const
{
  return &ivNegPatches;
}

void NNClassifier::removeWarps()
{
  for(unsigned int p = 0; p < ivPosPatches.size(); p++)
    if(ivWarpIndices[p] > 0)
    {
      ivPosPatches[p].erase(ivPosPatches[p].end() - ivWarpIndices[p], ivPosPatches[p].end());
      ivWarpIndices[p] = 0;
    }
}

void NNClassifier::addObject(const NNPatch& patch)
{
  ivPosPatches.push_back(std::vector<NNPatch>());
  ivPosPatches[ivPosPatches.size() - 1].push_back(patch);
  ivWarpIndices.push_back(0);
  // remove negative patches that are too similar to this new object
  for(int i = ivNegPatches.size() - 1; i >= 0; i--)
  {  
    double ncc = crossCorr(ivNegPatches[i].patch.data(), patch.patch.data(), 
                           ivNegPatches[i].norm2 * patch.norm2);
    if(ncc > 0.8){
      ivNegPatches.erase(ivNegPatches.begin() + i);    
      #if DEBUG
      std::cout << "removed negative patch " << i << " (ncc = " << ncc << ")" << std::endl;
      #endif
    }
  }
}

/// @param patch the patch that shall be learned
/// @param objId id of the object class to which the patch should be added
/// @param positive @b true if element of this class, @b false if not (background patch)
/// @param tmp @b true if this is a temporary (warped) patch
/// @returns @b true if the patch was added.
bool NNClassifier::trainNN(const NNPatch& patch, int objId, bool positive, bool tmp)
{
  double conf = getConf(patch, positive ? objId : -1, false);
  if(positive)
  {
    if(conf < 0.75)
    {
      ivPosPatches[objId].push_back(patch);    
      if(tmp)
        ivWarpIndices[objId]++;
      return true;
    }
  }
  else if(conf < 0.85)
  {
    ivNegPatches.push_back(patch);
    return true;
  }
  return false;
}

/// @param patch the patch that shall be evaluated
/// @param objId Id of the class to compare or -1 for comparison to negative (background) patches
/// @param conservative If @b true earlier positive patches are weighted more.
double NNClassifier::getConf(const NNPatch& patch, int objId, bool conservative) const
{
  if(ivUseColor)
  {
    double conf = getConf(patch.patch.data(), patch.norm2, objId, conservative);
    if(objId < 0 || conf < 0.58 || patch.histogram == NULL || ivPosPatches[objId][0].histogram == NULL)
      return conf;
    double colorCons = cmpHistograms(ivPosPatches[objId][0].histogram, patch.histogram);
        //Histogram::compareColorDistribution(ivPosPatches[objId][0].histogram, patch.histogram);
    //return std::min(colorCons, conf);
    if(colorCons < 0.75)
      return conf * 0.8;
    return conf + (1-conf) * 0.4;
  } //else
  return getConf(patch.patch.data(), patch.norm2, objId, conservative);
}

/// @see getConf(const NNPatch& patch, int objId, bool conservative)
double NNClassifier::getConf(NNPatch& patch, int objId, bool conservative,
                  const ObjectBox& bbox, const unsigned char * rgb, int w, int h) const
{
  if(ivUseColor)
  {
    double conf = getConf(patch.patch.data(), patch.norm2, objId, conservative);
    if(objId < 0 || conf < 0.58 || ivPosPatches[objId][0].histogram == NULL)
      return conf;
    else{
      if(patch.histogram == NULL && rgb != NULL)
        patch.histogram = Histogram::getInstance()->getColorDistribution(rgb, w, h, bbox);
      if(patch.histogram != NULL)
      {
        double colorCons = cmpHistograms(ivPosPatches[objId][0].histogram, patch.histogram);
          //Histogram::compareColorDistribution(ivPosPatches[objId][0].histogram, patch.histogram);
        //return std::min(colorCons, conf);
        if(colorCons < 0.75)
          return conf * 0.8;
        return conf + (1-conf) * 0.4;
      }
      return conf;
    }
  } //else
  return getConf(patch.patch.data(), patch.norm2, objId, conservative);
}

double NNClassifier::getConf(const float* patch, float norm2, int objId, bool conservative) const
{
  //max NCC with positive examples
  double posNCC = 0;
  if (objId >= 0)
  {
    int nPos = ivPosPatches[objId].size();
    double* posNCCs = new double[nPos];
    #pragma omp parallel for
    for (int i = 0; i < nPos; i++)
    {
      posNCCs[i] = crossCorr(ivPosPatches[objId][i].patch.data(), patch, (ivPosPatches[objId][i].norm2 * norm2));
      if (!ivAllowFastChange && conservative && i > nPos/2)
        posNCCs[i] *= 1.0 - 0.05 * (i-nPos/2) / (double)nPos;  
    }
    for (int i = 0; i < nPos; i++)
      if (posNCCs[i] > posNCC)
        posNCC = posNCCs[i]; 
    delete[] posNCCs;
  }
  //max NCC with negative examples
  double negNCC = 0;
  int nNeg = ivNegPatches.size();
  if (nNeg)
  {
    double* negNCCs = new double[nNeg];
    #pragma omp parallel for
    for (int i = 0; i < nNeg; i++)
      negNCCs[i] = crossCorr(ivNegPatches[i].patch.data(), patch, (ivNegPatches[i].norm2 * norm2));
    for (int i = 0; i < nNeg; i++)
      if (negNCCs[i] > negNCC)
        negNCC = negNCCs[i];
    delete[] negNCCs;
  }else
    negNCC = 0.3; //hack! there should always be a negative example from initialization
    
  if(objId < 0)
    return negNCC;
  return (1 - negNCC) / (2 - negNCC - posNCC);
}

double NNClassifier::crossCorr(const float* patchA, const float* patchB, float denom) const
{
  double sumDiff = 0;
  for (int i = 0; i < ivPatchSize*ivPatchSize; ++i)
    sumDiff += patchA[i] * patchB[i];
  if (denom <= 0)
    return sumDiff;
  return (sumDiff/sqrt(denom) + 1) / 2.0;
}

#define INVBINS (1.0/NUM_BINS)
double NNClassifier::cmpHistograms(const float* h1, const float* h2) const
{
  double corr = 0;
	double norm1 = 0;
	double norm2 = 0;
	for (int i = 0; i < NUM_BINS; ++i) 
	{
		corr  += (h1[i]-INVBINS) * (h2[i]-INVBINS);
		norm1 += (h1[i]-INVBINS) * (h1[i]-INVBINS);
		norm2 += (h2[i]-INVBINS) * (h2[i]-INVBINS);
	}
	return (corr / sqrt(norm1*norm2) + 1) / 2.0;
}

void NNClassifier::saveToStream(std::ofstream & outputStream) const
{
  outputStream.write((char*)&ivWidth, sizeof(int));
  outputStream.write((char*)&ivHeight, sizeof(int));
  outputStream.write((char*)&ivPatchSize, sizeof(int));
  outputStream.write((char*)&ivUseColor, sizeof(bool));
  outputStream.write((char*)&ivAllowFastChange, sizeof(bool));
  // negative patches
  int nNeg = ivNegPatches.size();
  outputStream.write((char*)&nNeg, sizeof(int));
  for(std::vector<NNPatch>::const_iterator it = ivNegPatches.begin(); it != ivNegPatches.end(); ++it)
    it->saveToStream(outputStream);
  // positive patches
  int nObs = ivPosPatches.size();
  outputStream.write((char*)&nObs, sizeof(int));
  for(std::vector<std::vector<NNPatch> >::const_iterator oit = ivPosPatches.begin(); oit != ivPosPatches.end(); ++oit)
  {
    int nPos = oit->size();
    outputStream.write((char*)&nPos, sizeof(int));
    for(std::vector<NNPatch>::const_iterator it = oit->begin(); it != oit->end(); ++it)
      it->saveToStream(outputStream);
  }
}

#endif //NNCLASSIFIER_H
