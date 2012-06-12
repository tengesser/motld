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


#ifndef FERNFILTER_H
#define FERNFILTER_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <vector>
#include <map>

#ifdef WIN32
#include <time.h>
#else
#include <sys/time.h>
#endif

#include "Matrix.h"
#include "Utils.h"

#define USEMAP 1       // Default: 1 - 0 = use lookup table instead: experimental
#define USETBBP 1      // Default: 1 - 0 = use simple pixel comparison: experimental
#define USEFASTSCAN 0  // Default: 0 - 1 = scan only every 2nd box: experimental

#define TIMING 1

// some settings - don't change these!
#define CONFIDENCETHRESHOLD             0.7
#define POSOVERLAPTHRESHOLD             0.85
#define NEGOVERLAPTHRESHOLD             0.65
#define INITNEGOVERLAPTHRESHOLD         0.20
#define NUMNEGTRAININGEXAMPLES         100
#define VARIANCETHRESHOLDFACTOR        0.8
#define VARIANCETHRESHOLDDESCENDRATE   0.2
#define VARIANCEMINTHRESHOLD           100


/// defines settings for affine warps, which are used in the FernFilter update process
struct WarpSettings
{
  /// (currently unused)
  int num_closest;
  /// how many warps should be done
  int num_warps;
  /// (currently unused)
  float noise;
  /// max rotation angle in degrees
  float angle;
  /// max translation shift in per cent
  float shift;
  /// max relative scale in per cent
  float scale;
};

/// describes a detection by the FernFilter
struct FernDetection
{
  /// the object box (location, dimension, object id)
  ObjectBox box;
  /// patch as downscaled, squared image
  Matrix patch;
  /// confidence value associated with patch
  float confidence;
  /// fern features associated with patch
  int * featureData;
  // temporary values
  const void * ss;     // pointer to scan parameters used for this detection
  float * imageOffset; // pointer to image / sat position required to compute featureData
  /// ordering over FernDetections (using their confidence values)
  static bool fdBetter(FernDetection fd1, FernDetection fd2) { return fd1.confidence > fd2.confidence; }
};

/// used for learning and (re-)finding objects
class FernFilter
{
public:
  /// public constructor
  FernFilter(const int & width, const int & height, const int & numFerns,
             const int & featuresPerFern, const int & patchSize = 15, 
             const int & scaleMin = -10, const int & scaleMax = 11, const int & bbMin = 24);
  /// copy constructor
  FernFilter(const FernFilter & other);
  /// destructor
  ~FernFilter();
  /// introduces new objects from a list of object boxes and returns negative training examples
  const std::vector<Matrix> addObjects(const Matrix & image, const std::vector<ObjectBox>& boxes);
  /// scans fern structure for possible object matches using a sliding window approach
  const std::vector<FernDetection> scanPatch(const Matrix & image) const;
  /// updates the fern structure with information about the correct boxes
  const std::vector< Matrix > learn(const Matrix& image, const std::vector< ObjectBox >& boxes, bool onlyVariance = false);
  /// creates a FernFilter from binary stream (load procedure)
  static FernFilter loadFromStream(std::ifstream & inputStream);
  /// writes FernFilter into binary stream (save procedure)
  void saveToStream(std::ofstream & outputStream) const;
  /// changes input image dimensions (has to be applied with applyPreferences())
  void changeInputFormat(const int & width, const int & height);
  /// changes default size scan box dimensions (has to be applied with applyPreferences())
  void changeScanBoxFormat(const int & width, const int & height);
  /// changes sliding window preferences (has to be applied with applyPreferences())
  void changeScanSettings(const int & scaleMin, const int & scaleMax, const int & bb_min);
  /// applies changes made by changeInputFormat() changeScanBoxFormat(), changeScanSettings()
  void applyPreferences();
  /// changes settings for warping
  void changeWarpSettings(const WarpSettings & initSettings, const WarpSettings & updateSettings);
  
private:
  // Methods for feature extraction / fern manipulation etc.
  void createScaledMatrix(const Matrix& image, Matrix & scaled, float*& sat, float*& sat2, int scale) const;
  void createScaledMatrices(const Matrix& image, Matrix*& scaled, float**& sats, float**& sat2s) const;
  void varianceFilter(float * image, float * sat, float * sat2, int scale, std::vector<FernDetection> & acc) const;
  std::vector< Matrix > retrieveHighVarianceSamples(const Matrix& image, const std::vector< ObjectBox >& boxes);
  int* extractFeatures(const float * const imageOrSAT, int ** offsets) const;
  void extractFeatures(FernDetection & det) const;
  float calcMaxConfidence(int * features) const;
  float * calcConfidences(int * features) const;
  void addPatch(const int & objId, const int * const featureData, const bool & pos);
  void addPatch(const Matrix& scaledImage, const int& objId, const bool& pos);
  void addPatchWithWarps(const Matrix & image, const ObjectBox & box, const WarpSettings & ws, 
                         std::vector<Matrix> & op, const bool & pos, const bool & notOnlyVar = true);
  void addWarpedPatches(const Matrix & image, const ObjectBox & box, const WarpSettings & ws,
                        std::vector<Matrix> & op, const bool & pos);
  void clearLastDetections() const;
  
  // Methods for initialization
  int *** createFeatures();
  void initializeFerns();
  void computeOffsets();
  int ** computeOffsets(int width);
  void addObjectToFerns();
  
  // Helper
  int calcTableSize() const;
  void debugOutput() const;
  FernDetection copyFernDetection(const FernDetection & fd) const;
  
  struct Posteriors
  {
    int n;
    int p;
    float posterior;
  };

  struct Confidences
  { 
    float maxConf;
    std::map<int, Posteriors> posteriors;
  };

  struct ScanSettings
  {
    int width;
    int height;
    float boxw;
    float boxh;
    float pixw;
    float pixh;
    int * varianceIndizes;
    int ** offsets;
  };
  
  #if DEBUG && TIMING
  struct ScanTime
  {
    int time0ScaledImages;
    int time1Variance;
    int time2GenFeatures;
    int time3CoarseFilter;
    int time4FineFilter;
    int time5GenPatches;
  };
#endif
  
  // changeable input image dimensions
  int ivWidth;
  int ivHeight;
  
  // hard classifier configuration
  const int ivNumFerns;
  const int ivFeaturesPerFern;
  const int ivPatchSize;
  
  // helper classifier constants
  const int ivPatchSizeMinusOne;
  const int ivPatchSizeSquared;
  
  // changable scan settings
  int ivOriginalWidth;
  int ivOriginalHeight;
  int ivScaleMin;
  int ivScaleMax;
  int ivBBmin;
  WarpSettings ivInitWarpSettings;
  WarpSettings ivUpdateWarpSettings;
  
  // Fern Data
  int *** ivFeatures;
#if USEMAP
  std::map<int, Confidences> * ivFernForest;
#else
  std::vector<int**> ivNtable;
  std::vector<int**> ivPtable;
  std::vector<float**> ivTable;
  float **ivMaxTable;
#endif
  
  // further instance variables
  int ivNumObjects;
  int ivScanNoZoom;
  float ivVarianceThreshold;
  int ** ivPatchSizeOffsets;
  std::vector<ScanSettings> ivScans;
  std::vector<float> ivMinVariances;
  mutable std::vector<FernDetection> ivLastDetections;
  
  // some default structures
  static const WarpSettings cDefaultInitWarpSettings;
  static const WarpSettings cDefaultUpdateWarpSettings;
};

const WarpSettings FernFilter::cDefaultInitWarpSettings = {10, 20, 5, 20, 0.2, 0.2};
const WarpSettings FernFilter::cDefaultUpdateWarpSettings = {10, 10, 5, 20, 0.2, 0.2};


/*****************************************************************************
*                         public accessible stuff                            *
******************************************************************************/

FernFilter::FernFilter(const int & width, const int & height, const int & numFerns,
                       const int & featuresPerFern, const int & patchSize, 
                       const int & scaleMin, const int & scaleMax, const int & bbMin)
                      : ivWidth(width), ivHeight(height), ivNumFerns(numFerns),
                        ivFeaturesPerFern(featuresPerFern), ivPatchSize(patchSize),
                        ivPatchSizeMinusOne(patchSize-1), 
                        ivPatchSizeSquared(patchSize*patchSize),
                        ivScaleMin(scaleMin), ivScaleMax(scaleMax), ivBBmin(bbMin),
                        ivInitWarpSettings(cDefaultInitWarpSettings),
                        ivUpdateWarpSettings(cDefaultUpdateWarpSettings),
                        ivFeatures(createFeatures()), ivNumObjects(0),
                        ivVarianceThreshold(255*255)
{
  initializeFerns();
}

const std::vector<Matrix> FernFilter::addObjects(const Matrix& image, const std::vector<ObjectBox>& boxes)
{
  std::vector<Matrix> result;
  std::vector<Matrix> posResult; // IDEA: positive warps could also be returned
  
  if (boxes.empty())
    return result;
  
  if (ivNumObjects == 0)
  {
    ivOriginalHeight = boxes[0].height;
    ivOriginalWidth = boxes[0].width;
    computeOffsets();
    //debugOutput();
  }
  
  for (unsigned int i = 0; i < boxes.size(); ++i)
  {
    if (boxes[i].objectId != ivNumObjects + (int)i)
      std::cerr << "ERROR WRONG OBJECT ENUMERATION!" << std::endl;
    ivMinVariances.push_back(100000);
    addObjectToFerns();
    addPatchWithWarps(image, boxes[i], ivInitWarpSettings, posResult, true);
  }
  
  if (ivNumObjects == 0)
    result = retrieveHighVarianceSamples(image, boxes);
  
  ivNumObjects += boxes.size();
  
  return result;
}

const std::vector<FernDetection> FernFilter::scanPatch(const Matrix & image) const
{ 
  // Pipeline structure
  std::vector<FernDetection> varianceFiltered;
  std::vector<FernDetection> fernFiltered1;
  std::vector<FernDetection> result;
  
  // scaled images, summed area tables
  Matrix* scaled; float** sats; float** sat2s;
  
  clearLastDetections();
  
  if (ivNumObjects == 0)
    return result;
  
#if DEBUG && TIMING
ScanTime st;
#endif
  
#pragma omp parallel
{

#if DEBUG && TIMING
  #pragma omp master
  st.time0ScaledImages = getTime();
#endif
  
  // Step 0 - Precalculate Scaled Images / Summed Area Tables
  createScaledMatrices(image, scaled, sats, sat2s); // HIER
  
#if DEBUG && TIMING
  #pragma omp master
  st.time1Variance = getTime();
#endif
  
  // STEP 1 - Scan, Filter by Variance
#pragma omp barrier
#pragma omp for schedule(dynamic)
  for (unsigned int i = 0; i < ivScans.size(); ++i)
    varianceFilter(scaled[i].data(), sats[i], sat2s[i], i, varianceFiltered);
  
#if DEBUG && TIMING
#pragma omp master
  st.time2GenFeatures = getTime();
#endif
  
  // STEP 2 - Calculate Feature Data
#pragma omp barrier
#pragma omp for
  for (unsigned int i = 0; i < varianceFiltered.size(); ++i)
    extractFeatures(varianceFiltered[i]);
  
#if DEBUG && TIMING
#pragma omp master
  st.time3CoarseFilter = getTime();
#endif
  
  // STEP 3 - Coarse filtering By Fern
#pragma omp barrier
  float confidenceThreshold = CONFIDENCETHRESHOLD * ivNumFerns;
#pragma omp for
  for (unsigned int i = 0; i < varianceFiltered.size(); ++i)
  {
    FernDetection det = varianceFiltered[i];
    det.confidence = calcMaxConfidence(det.featureData);
    if (det.confidence >= confidenceThreshold)
#pragma omp critical
      fernFiltered1.push_back(det);
    else
      delete[] det.featureData;
  }
  
#if DEBUG && TIMING
#pragma omp master
  st.time4FineFilter = getTime();
#endif
  
  // STEP 4 - Fine filtering By Fern
#pragma omp barrier
  if (ivNumObjects == 1)
  {
#pragma omp single
    result = fernFiltered1;
  } 
  else
  {
#pragma omp for
    for (unsigned int i = 0; i < fernFiltered1.size(); ++i)
    {
      FernDetection det = fernFiltered1[i];
      float * confidences = calcConfidences(det.featureData);
      int n = 0;
      for (int nObject = 0; nObject < ivNumObjects; ++nObject)
      {
        if (confidences[nObject] > confidenceThreshold)
        {
          FernDetection nDetection = copyFernDetection(det);
          nDetection.box.objectId = nObject;
#pragma omp critical
          result.push_back(nDetection);
          ++n;
        }
      }
      delete[] confidences;
      delete[] det.featureData;
    }
  }
  
#if DEBUG && TIMING
#pragma omp master
  st.time5GenPatches = getTime();
#endif
  
  // STEP 5 finally add patches
#pragma omp barrier
#pragma omp for
  for (unsigned int i = 0; i < result.size(); ++i)
  {
    result[i].patch.copyFromFloatArray(result[i].imageOffset,((ScanSettings*)(result[i].ss))->width,ivPatchSize,ivPatchSize);
  }
}
  

#if DEBUG && TIMING
  int timefinished = getTime();
#endif
  
#if DEBUG && TIMING
  // convert absolute time to relative time
  st.time0ScaledImages = st.time1Variance     - st.time0ScaledImages;
  st.time1Variance     = st.time2GenFeatures  - st.time1Variance;
  st.time2GenFeatures  = st.time3CoarseFilter - st.time2GenFeatures;
  st.time3CoarseFilter = st.time4FineFilter   - st.time3CoarseFilter;
  st.time4FineFilter   = st.time5GenPatches   - st.time4FineFilter;
  st.time5GenPatches   = timefinished         - st.time5GenPatches;
#endif
  
  delete[] scaled;
  for (unsigned int i = 0; i < ivScans.size(); ++i)
  {
    delete[] sats[i];
    delete[] sat2s[i];
  }
  delete[] sats;
  delete[] sat2s;
  
  ivLastDetections = result;
 
#if DEBUG  
  std::cout << "Patch Filterig Pipeline: " << varianceFiltered.size() << " >> " << fernFiltered1.size() 
            << " >> " << result.size();
#if TIMING
  std::cout << " | Time: " << st.time0ScaledImages << ", " << st.time1Variance << ", " << st.time2GenFeatures
            << ", " << st.time3CoarseFilter << ", " << st.time4FineFilter << ", " << st.time5GenPatches;
#endif
  std::cout << std::endl;
#endif
  
  return result;
}

// TODO: Multiprozessor
const std::vector<Matrix> FernFilter::learn(const Matrix & image, const std::vector<ObjectBox>& boxes, bool onlyVariance)
{
#if DEBUG
  int tStart = getTime();
#endif

  std::vector<Matrix> result;
  
  int del = 0;
  
  bool * valid = new bool[ivNumObjects];
  ObjectBox * bx = new ObjectBox[ivNumObjects];
  memset(valid, 0, ivNumObjects * sizeof(bool));
  
  // reinitialize varianceThreshold
  float lastVariance = ivVarianceThreshold;
  ivVarianceThreshold = 100000;
  
  // CATEGORIZE BOXES AND LEARN POSITIVE EXAMPLES
  for (std::vector<ObjectBox>::const_iterator bi = boxes.begin(); bi < boxes.end(); ++bi)
  {
    valid[bi->objectId] = true;
    bx[bi->objectId] = *bi;
    addPatchWithWarps(image, *bi, ivUpdateWarpSettings, result, true, !onlyVariance);
  }
  
  // calculate final variance value
  for (int nObj = 0; nObj < ivNumObjects; ++nObj)
  {
    if (!valid[nObj])
    {
      ivVarianceThreshold = MIN(ivVarianceThreshold, ivMinVariances[nObj]);
    }
  }
  float varThresholdSoll = MAX(VARIANCETHRESHOLDFACTOR*ivVarianceThreshold, VARIANCEMINTHRESHOLD);
  float diff = varThresholdSoll - lastVariance;
  ivVarianceThreshold = diff > 0 ? varThresholdSoll : lastVariance + diff * VARIANCETHRESHOLDDESCENDRATE;
  
  // UN-LEARN NEGATIVE EXAMPLES
  if (!onlyVariance)
  {
    for (std::vector<FernDetection>::iterator fd = ivLastDetections.begin();
  fd < ivLastDetections.end(); ++fd)
    {
      int objId = fd->box.objectId;
      if (valid[objId] && rectangleOverlap(fd->box, bx[objId]) < NEGOVERLAPTHRESHOLD)
      {
  addPatch(objId, fd->featureData, false);
  // negative patches don't have to be warped!
  // addPatch(image, fd->box, ivUpdateWarpSettings, false);
  ++del;
      }
    }
  }
  delete[] valid;
  delete[] bx;
  clearLastDetections();

#if DEBUG
  int tEnd = getTime();
  std::cout << "Fern Learner: +" << boxes.size() << " , -" << del << " | Variance: " << ivVarianceThreshold
            << " | Time: " << tEnd - tStart << std::endl;
#endif
  
  return result;
}

void FernFilter::saveToStream(std::ofstream & outputStream) const
{
  // 1. simple instance variables
  outputStream.write((char*)&ivWidth, sizeof(int));
  outputStream.write((char*)&ivHeight, sizeof(int));
  outputStream.write((char*)&ivNumObjects, sizeof(int));
  outputStream.write((char*)&ivNumFerns, sizeof(int));
  outputStream.write((char*)&ivFeaturesPerFern, sizeof(int));
  outputStream.write((char*)&ivPatchSize, sizeof(int));
  outputStream.write((char*)&ivScaleMin, sizeof(int));
  outputStream.write((char*)&ivScaleMax, sizeof(int));
  outputStream.write((char*)&ivBBmin, sizeof(int));
  outputStream.write((char*)&ivOriginalWidth, sizeof(int));
  outputStream.write((char*)&ivOriginalHeight, sizeof(int));
  outputStream.write((char*)&ivInitWarpSettings, sizeof(WarpSettings));
  outputStream.write((char*)&ivUpdateWarpSettings, sizeof(WarpSettings));
  outputStream.write((char*)&ivVarianceThreshold, sizeof(float));
  
  // 2. fern features
  for (int nFern = 0; nFern < ivNumFerns; ++nFern)
    for (int nFernFeature = 0; nFernFeature < ivFeaturesPerFern; ++nFernFeature)
      outputStream.write((char*)(ivFeatures[nFern][nFernFeature]), 4*sizeof(int));
  
  // 3. fern structures
#if USEMAP
  for (int nFern = 0; nFern < ivNumFerns; ++nFern)
  {
    int mapSize = ivFernForest[nFern].size();
    outputStream.write((char*)&mapSize, sizeof(int));
    for (std::map<int, Confidences>::iterator it = ivFernForest[nFern].begin(); it != ivFernForest[nFern].end(); ++it)
    {
      outputStream.write((char*)&(it->first), sizeof(int));
      outputStream.write((char*)&(it->second.maxConf), sizeof(float));
      int entrySize = it->second.posteriors.size();
      outputStream.write((char*)&entrySize, sizeof(int));
      for (std::map<int, Posteriors>::iterator it2  = it->second.posteriors.begin(); 
                                         it2 != it->second.posteriors.end(); ++it2)
      {
        outputStream.write((char*)&(it2->first), sizeof(int));
        outputStream.write((char*)&(it2->second), sizeof(Posteriors));
      }
    }
  }
#else
  std::cerr << "Saving Not Yet implemented for Table" << std::endl;
#endif

  // 4. minVariance
  outputStream.write((char*)ivMinVariances.data(), ivNumObjects*sizeof(float));
    
  /* DEBUGGING - print out instance variables
  std::cout << ivWidth << ", " << ivHeight << ", " << ivNumObjects << ", " << ivNumFerns
            << ", " << ivFeaturesPerFern << ", " << ivPatchSize << ", " << ivScaleMin
            << ", " << ivScaleMax << ", " << ivBBmin << ", " << ivVarianceThreshold << std::endl;
  std::cout << ivInitWarpSettings.angle << ", " << ivInitWarpSettings.noise << ", " 
            << ivInitWarpSettings.scale << ", " << ivInitWarpSettings.shift << " | "
      << ivInitWarpSettings.num_closest << ", " << ivInitWarpSettings.num_warps << std::endl;
  std::cout << ivUpdateWarpSettings.angle << ", " << ivUpdateWarpSettings.noise << ", " 
            << ivUpdateWarpSettings.scale << ", " << ivUpdateWarpSettings.shift << " | "
      << ivUpdateWarpSettings.num_closest << ", " << ivUpdateWarpSettings.num_warps << std::endl;
   */
}

FernFilter FernFilter::loadFromStream(std::ifstream & inputStream)
{
  // 1. simple instance variables
  int width, height, numObjects, numFerns, featuresPerFern, patchSize, scaleMin, scaleMax,
      bbMin, originalWidth, originalHeight;
  WarpSettings initWarpSettings, updateWarpSettings;
  float varianceThreshold;
  inputStream.read((char*)&width, sizeof(int));
  inputStream.read((char*)&height, sizeof(int));
  inputStream.read((char*)&numObjects, sizeof(int));
  inputStream.read((char*)&numFerns, sizeof(int));
  inputStream.read((char*)&featuresPerFern, sizeof(int));
  inputStream.read((char*)&patchSize, sizeof(int));
  inputStream.read((char*)&scaleMin, sizeof(int));
  inputStream.read((char*)&scaleMax, sizeof(int));
  inputStream.read((char*)&bbMin, sizeof(int));
  inputStream.read((char*)&originalWidth, sizeof(int));
  inputStream.read((char*)&originalHeight, sizeof(int));
  inputStream.read((char*)&initWarpSettings, sizeof(WarpSettings));
  inputStream.read((char*)&updateWarpSettings, sizeof(WarpSettings));
  inputStream.read((char*)&varianceThreshold, sizeof(float));
  
  // 2. fern features
  int *** features = new int**[numFerns];
  for (int nFern = 0; nFern < numFerns; ++nFern)
  {
    features[nFern] = new int*[featuresPerFern];
    for (int nFernFeature = 0; nFernFeature < featuresPerFern; ++nFernFeature)
    {
      features[nFern][nFernFeature] = new int[4];
      inputStream.read((char*)(features[nFern][nFernFeature]), 4*sizeof(int));
    }
  }
  
  // 3. fern structures
#if USEMAP
  std::map<int, Confidences> * fernForest = new std::map<int, Confidences>[numFerns];
  for (int nFern = 0; nFern < numFerns; ++nFern)
  {
    int mapSize;
    inputStream.read((char*)&mapSize, sizeof(int));
    for (int nFernLeaf = 0; nFernLeaf < mapSize; ++nFernLeaf)
    {
      int key, entrySize; float maxConf;
      inputStream.read((char*)&key, sizeof(int));
      inputStream.read((char*)&maxConf, sizeof(float));
      inputStream.read((char*)&entrySize, sizeof(int));
      std::map<int, Posteriors> posteriors;
      for (int nEntry = 0; nEntry < entrySize; ++nEntry)
      {
  int nObject; Posteriors posterior;
  inputStream.read((char*)&nObject, sizeof(int));
        inputStream.read((char*)&posterior, sizeof(Posteriors));
        posteriors[nObject] = posterior;
      }
      Confidences cf;
      cf.maxConf = maxConf;
      cf.posteriors = posteriors;
      fernForest[nFern][key] = cf;
    }
  }
#else
  std::cerr << "Loading Not Yet implemented for Table" << std::endl;
#endif

  // 4. minVariance
  std::vector<float> minVariances(numObjects);
  inputStream.read((char*)minVariances.data(), numObjects * sizeof(float));
  
  // finally generate fern filter
  FernFilter result(width, height, numFerns, featuresPerFern, patchSize, scaleMin, scaleMax, bbMin);
  result.ivNumObjects = numObjects;
  result.ivInitWarpSettings = initWarpSettings;
  result.ivUpdateWarpSettings = updateWarpSettings;
  result.ivVarianceThreshold = varianceThreshold;
  result.ivOriginalWidth = originalWidth;
  result.ivOriginalHeight = originalHeight;
  result.ivFeatures = features;
#if USEMAP
  result.ivFernForest = fernForest;
#endif
  result.ivMinVariances = minVariances;
  result.computeOffsets();
  // result.debugOutput();
  
  /* DEBUGGING - print out loaded instance variables
  std::cout << width << ", " << height << ", " << numObjects << ", " << numFerns
            << ", " << featuresPerFern << ", " << patchSize << ", " << scaleMin
            << ", " << scaleMax << ", " << bbMin << ", " << varianceThreshold << std::endl;
  std::cout << initWarpSettings.angle << ", " << initWarpSettings.noise << ", " 
            << initWarpSettings.scale << ", " << initWarpSettings.shift << " | "
            << initWarpSettings.num_closest << ", " << initWarpSettings.num_warps << std::endl;
  std::cout << updateWarpSettings.angle << ", " << updateWarpSettings.noise << ", " 
            << updateWarpSettings.scale << ", " << updateWarpSettings.shift << " | "
            << updateWarpSettings.num_closest << ", " << updateWarpSettings.num_warps << std::endl;
  std::cout << result.ivFernForest->size() << ", " << result.ivScans.size() << std::endl;
  */
 
  return result;
}

FernFilter::FernFilter(const FernFilter& source) :
  ivWidth(source.ivWidth), ivHeight(source.ivHeight), 
  ivNumFerns(source.ivNumFerns), ivFeaturesPerFern(source.ivFeaturesPerFern),
  ivPatchSize(source.ivPatchSize), ivPatchSizeMinusOne(source.ivPatchSizeMinusOne), 
  ivPatchSizeSquared(source.ivPatchSizeSquared),
  ivOriginalWidth(source.ivOriginalWidth),
  ivOriginalHeight(source.ivOriginalHeight),
  ivScaleMin(source.ivScaleMin), 
  ivScaleMax(source.ivScaleMax), ivBBmin(source.ivBBmin),
  ivInitWarpSettings(source.ivInitWarpSettings),
  ivUpdateWarpSettings(source.ivUpdateWarpSettings),
  ivNumObjects(source.ivNumObjects),
  ivScanNoZoom(source.ivScanNoZoom),
  ivVarianceThreshold(source.ivVarianceThreshold),
  ivMinVariances(source.ivMinVariances)
{
  // copy ivFeatures
  ivFeatures = new int**[ivNumFerns];
  for (int nFern = 0; nFern < ivNumFerns; ++nFern)
  {
    ivFeatures[nFern] = new int*[ivFeaturesPerFern];
    for (int nFeature = 0; nFeature < ivFeaturesPerFern; ++nFeature)
    {
      ivFeatures[nFern][nFeature] = new int[4];
      memcpy(ivFeatures[nFern][nFeature], source.ivFeatures[nFern][nFeature], 4*sizeof(int));
    }
  }
  
  // copy fern data
#if USEMAP
  ivFernForest = new std::map<int, Confidences>[ivNumFerns];
  for (int nFern = 0; nFern < ivNumFerns; ++nFern)
  {
    ivFernForest[nFern] = source.ivFernForest[nFern];
  }
#else
  std::cerr << "COPY CONSTRUCTOR NOT YET IMPLEMENTED FOR LOOKUP TABLE!" << std::endl;
#endif
    
#if USEMAP
  int offsetSize = 16*ivFeaturesPerFern;
#else
  int offsetSize = 4*ivFeaturesPerFern;
#endif
    
  // copy scan data
  if (source.ivNumObjects > 0)
  {
    for (std::vector<ScanSettings>::const_iterator it = source.ivScans.begin(); 
                                             it < source.ivScans.end();
                                                   ++it)
    {
      ScanSettings ss = *it;
      ss.offsets = new int*[ivNumFerns];
      for (int nFern = 0; nFern < ivNumFerns; ++nFern)
      {
  ss.offsets[nFern] = new int[offsetSize];
  memcpy(ss.offsets[nFern], it->offsets[nFern], offsetSize * sizeof(int));
      }
      ss.varianceIndizes = new int[4];
      memcpy(ss.varianceIndizes, it->varianceIndizes, 4 * sizeof(int));
      ivScans.push_back(ss);
    }
    
    //ivPatchSizeOffsets;
    ivPatchSizeOffsets = new int*[ivNumFerns];
    for (int nFern = 0; nFern < ivNumFerns; ++nFern)
    {
      ivPatchSizeOffsets[nFern] = new int[offsetSize];
      memcpy(ivPatchSizeOffsets[nFern], source.ivPatchSizeOffsets[nFern], offsetSize * sizeof(int));
    }
  }
}

FernFilter::~FernFilter()
{
  // Learned Data
#if USEMAP
  delete [] ivFernForest;
#else
  std::cerr << "Destructor not implemented for table lookup!" << std::endl;
#endif
  
  if (ivNumObjects > 0)
  {
    // Scan Settings / Offsets
    for (std::vector<ScanSettings>::iterator it = ivScans.begin(); it < ivScans.end(); ++it)
    {
      delete[] it->varianceIndizes;
      for (int nFern = 0; nFern < ivNumFerns; ++nFern)
  delete[] it->offsets[nFern];
      delete[] it->offsets;
    }
    // Patch Size Offsets
    for (int nFern = 0; nFern < ivNumFerns; ++nFern)
    {
      delete[] ivPatchSizeOffsets[nFern];
    }
    delete[] ivPatchSizeOffsets;
  }
  // Features
  for (int nFern = 0; nFern < ivNumFerns; ++nFern)
  {
    for (int nFeature = 0; nFeature < ivFeaturesPerFern; ++nFeature)
    {
      delete[] ivFeatures[nFern][nFeature];
    }
    delete[] ivFeatures[nFern];
  }
  delete[] ivFeatures;
  clearLastDetections();
}

/*****************************************************************************
*                                   setters                                  *
******************************************************************************/

void FernFilter::changeInputFormat(const int& width, const int& height)
{
  ivWidth = width;
  ivHeight = height;
}

void FernFilter::changeScanBoxFormat(const int& width, const int& height)
{
  ivOriginalWidth = width;
  ivOriginalHeight = height;
}

void FernFilter::changeScanSettings(const int& scaleMin, const int& scaleMax, const int& bb_min)
{
  ivScaleMin = scaleMin;
  ivScaleMax = scaleMax;
  ivBBmin = bb_min;
}

void FernFilter::applyPreferences()
{
  std::cout << "Applying Input and Scanning settings..." << std::endl;
  computeOffsets();
  // debugOutput();
}

void FernFilter::changeWarpSettings(const WarpSettings& initSettings, const WarpSettings& updateSettings)
{
  ivInitWarpSettings = initSettings;
  ivUpdateWarpSettings = updateSettings;
}

/*****************************************************************************
*                         private accessible stuff                           *
******************************************************************************/

inline void FernFilter::createScaledMatrix(const Matrix& image, Matrix& scaled, float*& sat, float*& sat2, int scale) const
{
  ScanSettings ss = ivScans[scale];
  scaled = image;
  scaled.rescale(ss.width, ss.height);
  float ** saTables = scaled.createSummedAreaTable2();
  sat  = saTables[0];
  sat2 = saTables[1];
  delete[] saTables;
}

inline void FernFilter::createScaledMatrices(const Matrix& image, Matrix*& scaled, float**& sats, float**& sat2s) const
{
#pragma omp master
{
  const int numScans = ivScans.size();
  scaled = new Matrix[numScans];
  sats   = new float*[numScans];
  sat2s  = new float*[numScans];
}
#pragma omp barrier
#pragma omp for schedule(dynamic)
  for (unsigned int i = 0; i < ivScans.size(); ++i)
  {
    createScaledMatrix(image, scaled[i], sats[i], sat2s[i], i);
  }
}

inline void FernFilter::varianceFilter(float * image, float * sat, float * sat2, int scale, std::vector<FernDetection> & acc) const
{   
  ScanSettings ss = ivScans[scale];
  int right = ss.width - ivPatchSizeMinusOne;
  int bottom = ss.height - ivPatchSizeMinusOne;

  for (int y = 0; y < bottom; ++y)
  {
    int yDiff = y * (ss.width + 1);
    float * satPos = sat + yDiff;
    float * sat2Pos = sat2 + yDiff;
    float * imgPos = image + y * ss.width;
      
#if USEFASTSCAN
    int fst = y % 2;
    if (fst == 1) { imgPos++; satPos++; sat2Pos++; }
    int step = 2;
#else
    int fst = 0;
    int step = 1;   
#endif
      
    for (int x = fst; x < right; x += step, sat2Pos += step, satPos += step, imgPos += step) //, ++start
    {
      float ex2 = summedTableArea(sat2Pos,ss.varianceIndizes)/ivPatchSizeSquared;
      float ex = summedTableArea(satPos,ss.varianceIndizes)/ivPatchSizeSquared;
      float variance = ex2 - ex*ex;
        
      if (variance >= ivVarianceThreshold)
      {
        FernDetection fd;
        ObjectBox box = {x*ss.pixw, y*ss.pixh, ss.boxw, ss.boxh, 0};
        fd.imageOffset = imgPos;
        fd.box = box;
        fd.confidence = variance;
        fd.ss = &(ivScans[scale]);
        
#if USETBBP
        fd.featureData = (int*)satPos;
#else
        fd.featureData = (int*)imgPos;
#endif
#pragma omp critical
        acc.push_back(fd);
      }
    }
  }
}

inline std::vector<Matrix> FernFilter::retrieveHighVarianceSamples(const Matrix& image, const std::vector<ObjectBox>& boxes)
{
  std::vector<Matrix> result;
  
  // generate Summed Area Tables
  Matrix scaled; 
  float* sat;
  float* sat2;
  createScaledMatrix(image, scaled, sat, sat2, ivScanNoZoom);
  
  // scan and order hits
  std::vector<FernDetection> varianceDetections;
  float ivVarTTmp = ivVarianceThreshold;
  ivVarianceThreshold = VARIANCEMINTHRESHOLD;
  varianceFilter(scaled.data(), sat, sat2, ivScanNoZoom, varianceDetections);
  ivVarianceThreshold = ivVarTTmp;
  std::sort(varianceDetections.begin(), varianceDetections.end(), FernDetection::fdBetter);
  
  std::vector<ObjectBox> allBoxes = boxes;
  for (unsigned int i = 0; i < varianceDetections.size() && result.size() < NUMNEGTRAININGEXAMPLES; ++i)
  {
    bool good = true;
    for (unsigned int j = 0; good && j < allBoxes.size(); ++j)
      if (rectangleOverlap(varianceDetections[i].box, allBoxes[j]) > INITNEGOVERLAPTHRESHOLD)
        good = false;
    if (good)
    {
      allBoxes.push_back(varianceDetections[i].box);
      Matrix m;
      m.copyFromFloatArray(varianceDetections[i].imageOffset, ivScans[ivScanNoZoom].width, ivPatchSize, ivPatchSize);
      result.push_back(m);
    }
  }
  
  delete[] sat;
  delete[] sat2;
  
  return result;
}

inline int*** FernFilter::createFeatures()
{
  int *** result = new int**[ivNumFerns];
  for (int nFern = 0; nFern < ivNumFerns; ++nFern)
  {
    result[nFern] = new int*[ivFeaturesPerFern];
    for (int nFeature = 0; nFeature < ivFeaturesPerFern; ++nFeature)
    {
      result[nFern][nFeature] = new int[4];
      result[nFern][nFeature][2] = randInt(2, ivPatchSize); // width
      result[nFern][nFeature][3] = randInt(2, ivPatchSize); // height
      result[nFern][nFeature][0] = randInt(0, ivPatchSize - result[nFern][nFeature][2]); // x position
      result[nFern][nFeature][1] = randInt(0, ivPatchSize - result[nFern][nFeature][3]); // y position
    }
  }
  return result;
}

inline void FernFilter::initializeFerns()
{
#if USEMAP
  ivFernForest = new std::map<int, Confidences>[ivNumFerns];
  
#else
  int tableSize = calcTableSize();
  
  ivMaxTable = new float*[ivNumFerns];
  for (int nFern = 0; nFern < ivNumFerns; ++nFern)
  {
    ivMaxTable[nFern] = new float[tableSize];
    memset(ivMaxTable[nFern], 0, tableSize * sizeof(float));
  }
#endif
}

inline void FernFilter::addObjectToFerns()
{
#if !USEMAP
  int tableSize = calcTableSize();
  int ** nTable = new int*[ivNumFerns];
  int ** pTable = new int*[ivNumFerns];
  float ** table = new float*[ivNumFerns];
  for (int nFern = 0; nFern < ivNumFerns; ++nFern)
  {
    nTable[nFern] = new int[tableSize];
    memset(nTable[nFern], 0, tableSize*sizeof(int));
    pTable[nFern] = new int[tableSize];
    memset(pTable[nFern], 0, tableSize*sizeof(int));
    table[nFern] = new float[tableSize];
    memset(table[nFern], 0., tableSize*sizeof(int));
  }
  ivNtable.push_back(nTable);
  ivPtable.push_back(pTable);
  ivTable.push_back(table);
#endif
}

inline void FernFilter::computeOffsets()
{
  ivScans.clear();
  ivScanNoZoom  = 0;
  
  for (int scli = ivScaleMin; scli < ivScaleMax; ++scli)
  {
    ScanSettings ss;
    
    float scale = pow(1.2, scli);
    
    ss.boxw = ivOriginalWidth * scale;
    ss.boxh = ivOriginalHeight * scale;
    
    if (ss.boxw < ivBBmin || ss.boxh < ivBBmin || ss.boxw > ivWidth || ss.boxh > ivHeight)
      continue;
    
    ss.pixw = ss.boxw / ivPatchSize;
    ss.pixh = ss.boxh / ivPatchSize;
    
    ss.width  = round(ivWidth  / ss.pixw);
    ss.height = round(ivHeight / ss.pixh);
    
    ss.varianceIndizes = getSATIndices(ss.width, ivPatchSize, ivPatchSize);
    ss.offsets = computeOffsets(ss.width);
    
    if (scli == 0)
      ivScanNoZoom = ivScans.size();
    
    ivScans.push_back(ss);
  }
  
  ivPatchSizeOffsets = computeOffsets(ivPatchSize);
}

inline int ** FernFilter::computeOffsets(int width)
{
#if USETBBP
  int ** result = new int*[ivNumFerns];
  for (int nFern = 0; nFern < ivNumFerns; ++nFern)
  {
    int * fernOffsets = new int[ivFeaturesPerFern*16];
    for (int nFeature = 0; nFeature < ivFeaturesPerFern; ++nFeature)
    {
      int * feature = ivFeatures[nFern][nFeature];
      int * currentPos = fernOffsets + 16 * nFeature;
      int hpix = feature[2] >> 1;
      int left = feature[0]; // + x = + 0
      int right = left + feature[2] - 1;
      int centerl = left + hpix - 1;
      int centerr = right - hpix + 1;
    
      int vpix = feature[3] >> 1;
      int top = feature[1]; // + y = + 0
      int bottom = top + feature[3] - 1;
      int middlet = top + vpix - 1;
      int middleb = bottom - vpix + 1;
      
      getSATIndices(currentPos     , width, left   ,     top,   right, middlet);
      getSATIndices(currentPos +  4, width, left   , middleb,   right,  bottom);
      getSATIndices(currentPos +  8, width, left   ,     top, centerl,  bottom);
      getSATIndices(currentPos + 12, width, centerr,     top,   right,  bottom);
    }
    result[nFern] = fernOffsets;
  }
  return result;
#else
  int ** result = new int*[ivNumFerns];
  for (int nFern = 0; nFern < ivNumFerns; ++nFern)
  {
    int * fernOffsets = new int[ivFeaturesPerFern*2];
    int i = 0;
    for (int nFeature = 0; nFeature < ivFeaturesPerFern; ++nFeature)
    {
      int* feature = ivFeatures[nFern][nFeature];
      int p1 = feature[1] * width + feature[0];
      int p2 = p1 + feature[3] * width + feature[2];
      fernOffsets[i++] = p1;
      fernOffsets[i++] = p2;
    }
    result[nFern] = fernOffsets;
  }
  return result;
#endif
}

inline int FernFilter::calcTableSize() const
{
  int result = 1 << ivFeaturesPerFern;
 #if USETBBP
  result *= result;
 #endif
  return result;
}

inline void FernFilter::debugOutput() const
{
  // FEATURES
  for (int nFern = 0; nFern < ivNumFerns; ++nFern)
  {
    std::cout << "Fern " << nFern << ":" << std::endl;
    for (int i = -1; i < 4; ++i)
    {
      std::cout << (i == - 1 ? "feat" : i == 0 ? "posX" : i == 1 ? "posY" : i == 2 ? "lenX" : "lenY") << "\t";
      for (int nFeature = 0; nFeature < ivFeaturesPerFern; ++nFeature)
      {
        std::cout << (i == -1 ? nFeature : ivFeatures[nFern][nFeature][i]) << "\t";
      }
      std::cout << std::endl;
    }
  }
 
  // SCAN SETTINGS
  std::cout << "Scan Settings\nboxw\tboxh\twdth\thght\tpixw\tpixh\tvi1\tvi2\tvi3\tvi4" << std::endl;
  
  for (std::vector<ScanSettings>::const_iterator it = ivScans.begin(); it < ivScans.end(); ++it)
  {
    std::cout << it->boxw << "\t" << it->boxh << "\t" << it->width << "\t" << it->height << "\t" 
              << it->pixw << "\t" << it->pixh << "\t" << it->varianceIndizes[0] << "\t" <<
                 it->varianceIndizes[1] << "\t" << it->varianceIndizes[2] << "\t" <<
                 it->varianceIndizes[3] << std::endl;
  }
  
}

inline float FernFilter::calcMaxConfidence(int* featureData) const
{
  float result = 0;
  for (int nFern = 0; nFern < ivNumFerns; ++nFern)
  {
#if USEMAP
      std::map<int, Confidences>::const_iterator found =
        ivFernForest[nFern].find(featureData[nFern]);
      if (found != ivFernForest[nFern].end())
      {
        result += found->second.maxConf;
      }
#else
      int f = featureData[nFern];
      result += ivMaxTable[nFern][f];
#endif
  }
  return result;
}

inline float * FernFilter::calcConfidences(int* features) const
{
  float * result = new float[ivNumObjects];
  memset(result, 0, ivNumObjects * sizeof(float));
#if USEMAP
  for (int nFern = 0; nFern < ivNumFerns; ++nFern)
  {
    std::map<int, Confidences>::const_iterator found = 
      ivFernForest[nFern].find(features[nFern]);
    if (found != ivFernForest[nFern].end())
    {
      for (std::map<int, Posteriors>::const_iterator pi = found->second.posteriors.begin();
           pi != found->second.posteriors.end(); ++pi)
      {
        result[pi->first] += pi->second.posterior;
      }
    }
  }
#else
  for (int nObject = 0; nObject < ivNumObjects; ++nObject)
  {
    float ** table = ivTable[nObject];
    for (int nFern = 0; nFern < ivNumFerns; ++nFern)
    {
      result[nObject] += table[nFern][features[nFern]];
    }
  }
#endif
  return result;
}

inline void FernFilter::extractFeatures(FernDetection& det) const
{
  det.featureData = extractFeatures((float*)det.featureData, ((ScanSettings*)(det.ss))->offsets);
}

inline int * FernFilter::extractFeatures(const float * const imgOrSAT, int ** offsets) const
{
  int * result = new int[ivNumFerns];
  for(int nFern = 0; nFern < ivNumFerns; ++nFern)
  {
    int * foffsets = offsets[nFern];
#if USETBBP
    int fernClass = 0;
    for (int nFeature = 0; nFeature < ivFeaturesPerFern; ++nFeature)
    {
      int * foffsets2 = foffsets + 16 * nFeature;
      float artop = summedTableArea(imgOrSAT, foffsets2);
      float arbot = summedTableArea(imgOrSAT, foffsets2+4);
      float arlft = summedTableArea(imgOrSAT, foffsets2+8);
      float arrgt = summedTableArea(imgOrSAT, foffsets2+12);
      int vf = artop < arbot;
      int hf = arlft < arrgt;
      fernClass <<= 2;
      fernClass |= (vf << 1) | hf;
    }
#else
    int fernClass = 0;
    for (int nFeature = 0; nFeature < ivFeaturesPerFern; ++nFeature)
    {
      float px1 = imgOrSAT[*foffsets++];
      float px2 = imgOrSAT[*foffsets++];
      fernClass = (fernClass << 1) | (px1 < px2);
    }
#endif
    result[nFern] = fernClass;
  }
  return result;
}

inline void FernFilter::addPatch(const int & objId, const int * const featureData, const bool & pos)
{
#if !USEMAP
  int ** nTable = ivNtable[objId];
  int ** pTable = ivPtable[objId];
  float ** table = ivTable[objId];
#endif
  for (int nFern = 0; nFern < ivNumFerns; ++nFern)
  {
    // TODO: nFern ausgliedern?
    int feature = featureData[nFern];
#if USEMAP
    std::map<int, Confidences>::iterator found =
        ivFernForest[nFern].find(feature);
    if (found != ivFernForest[nFern].end())
    {
      std::map<int, Posteriors>::iterator found2 =
        found->second.posteriors.find(objId);
      if (found2 != found->second.posteriors.end())
      {
        bool isPos = found2->second.p >= CONFIDENCETHRESHOLD;
        if (pos != isPos)
        {
          (pos ? found2->second.p : found2->second.n) += 1;
          int pf = found2->second.p;
          found2->second.posterior = (float)pf / (pf + found2->second.n);
          found->second.maxConf = 0;
          for (std::map<int, Posteriors>::const_iterator it = found->second.posteriors.begin();
               it != found->second.posteriors.end(); ++it)
          {
            found->second.maxConf = MAX(found->second.maxConf, it->second.posterior);
          }
        }
      }
      else
      {
        Posteriors p = {pos ? 0 : 1, pos ? 1 : 0, pos ? 1 : 0};
        found->second.posteriors[objId] = p;
        found->second.maxConf = MAX(found->second.maxConf, p.posterior);
      }
    }
    else
    {
      Posteriors p = {pos ? 0 : 1, pos ? 1 : 0, pos ? 1 : 0};
      Confidences c;
      c.posteriors[objId] = p;
      c.maxConf = p.posterior;
      ivFernForest[nFern][feature] = c;
    }
#else
    bool isPos = table[nFern][feature] >= CONFIDENCETHRESHOLD;
    if (pos != isPos)
    {
      (pos ? pTable : nTable)[nFern][feature]++;
      int pf = pTable[nFern][feature];
      table[nFern][feature] = (float)pf / (pf + nTable[nFern][feature]);
      ivMaxTable[nFern][feature] = 0;
      for (int numObj = 0; numObj < ivNumObjects; ++numObj)
      {
        ivMaxTable[nFern][feature] = MAX(ivMaxTable[nFern][feature], ivTable[numObj][nFern][feature]);
      }
    }
#endif
  }
}

inline void FernFilter::addPatch(const Matrix& scaledImage, const int& objId, const bool& pos)
{
#if USETBBP
  float * sat = scaledImage.createSummedAreaTable();
  int * features = extractFeatures(sat, ivPatchSizeOffsets);
  delete[] sat;
#else
  int * features = extractFeatures(scaledImage.data(), ivPatchSizeOffsets);
#endif
  addPatch(objId, features, pos);
  delete[] features;
}

inline void FernFilter::addPatchWithWarps(const Matrix& image, const ObjectBox& box, const WarpSettings& ws, 
                                          std::vector<Matrix> & op, const bool& pos, const bool& notOnlyVar)
{ 
  float factX = box.width / ivPatchSize;
  float factY = box.height / ivPatchSize;
  float width = image.xSize() / factX;
  float height = image.ySize() / factY;
  
  Matrix scaled = image; scaled.rescale(round(width), round(height));  
  ObjectBox newB = {box.x / factX, box.y / factY, ivPatchSize, ivPatchSize};
  newB.objectId = box.objectId;
  
  Matrix pt(box.width,box.height);
  pt.copyFromFloatArray(scaled.data(), width, height, round(newB.x), round(newB.y), round(newB.width), round(newB.height));
  pt.rescale(ivPatchSize,ivPatchSize);
  
  float ** sats = pt.createSummedAreaTable2();
  const int index = (ivPatchSize+1)*(ivPatchSize+1)-1;
  const int nPixels = ivPatchSize * ivPatchSize;
  const float ex2 = sats[1][index] / nPixels;
  const float ex = sats[0][index] / nPixels;
  const float exex = ex * ex;
  const float variance = ex2 - exex;
  
  ivMinVariances[box.objectId] = MIN(ivMinVariances[box.objectId], variance);
  ivVarianceThreshold = MIN(ivVarianceThreshold, variance);
  
  if (notOnlyVar)
  {
#if USETBBP
    float * img = sats[0];
#else
    float * img = pt.data();
#endif
    const int * data = extractFeatures(img, ivPatchSizeOffsets);
    addPatch(box.objectId, data, pos);
    addWarpedPatches(scaled, newB, ws, op, pos);
    delete[] data;
  }
  
  delete[] sats[0];
  delete[] sats[1];
  delete[] sats;
}

inline void FernFilter::addWarpedPatches(const Matrix& image, const ObjectBox& box, const WarpSettings& ws, 
                                         std::vector<Matrix> & op, const bool& pos)
{  
  
  // Default Warps!
  Matrix defWl = Matrix::createWarpMatrix( 0.2, 1);
  //Matrix defWn = createWarpMatrix(   0, 1);
  Matrix defWr = Matrix::createWarpMatrix(-0.2, 1);
  
  Matrix imgWarpL = image.affineWarp(defWl, box, false);
  //Matrix imgWarpN = image.affineWarp(defWn, box, false);
  Matrix imgWarpR = image.affineWarp(defWr, box, false);
  
  /* DEBUGGING - print affine warps
  static int x = 0;
  char filenameL[255];
  char filenameN[255];
  char filenameR[255];
  sprintf(filenameL, "output/warped_pic%03d_obj%d_l.pgm",x, box.objectId);
  sprintf(filenameN, "output/warped_pic%03d_obj%d_n.pgm",x, box.objectId);
  sprintf(filenameR, "output/warped_pic%03d_obj%d_r.pgm",x, box.objectId);
  imgWarpL.writeToPGM(filenameL);
  imgWarpN.writeToPGM(filenameN);
  imgWarpR.writeToPGM(filenameR);
  ++x;
  */
  
  
  imgWarpL.rescale(ivPatchSize, ivPatchSize);
  imgWarpR.rescale(ivPatchSize, ivPatchSize);
  addPatch(imgWarpL, box.objectId, pos);
  addPatch(imgWarpR, box.objectId, pos);
  op.push_back(imgWarpL);
  op.push_back(imgWarpR);
  
  for (int i = 0; i < ws.num_warps; ++i)
  {
    float angle = (PI / 180) * ws.angle * randFloat(-0.5, 0.5);
    float scale = 1 - ws.scale * randFloat(-0.5, 0.5);  
    Matrix warpMatrix = Matrix::createWarpMatrix(angle, scale);
    
    float shX = ws.shift * box.width * randFloat(-0.5, 0.5);
    float shY = ws.shift * box.height * randFloat(-0.5, 0.5);
    ObjectBox shiftedBox = {box.x + shX, box.y + shY, box.width, box.height};
    
    Matrix warped = image.affineWarp(warpMatrix, shiftedBox, true);
    
    /*char filename[255];
    sprintf(filename, "output/warped2_pic%03d_obj%d_%d.pgm", x, box.objectId, i);
    warped.writeToPGM(filename);*/
    
    warped.rescale(ivPatchSize,ivPatchSize);
    addPatch(warped, box.objectId, pos);
  }
}

inline void FernFilter::clearLastDetections() const
{
  for (std::vector<FernDetection>::iterator it = ivLastDetections.begin(); it < ivLastDetections.end(); ++it)
  {
    delete[] it->featureData;
  }
  ivLastDetections.clear();
}

inline FernDetection FernFilter::copyFernDetection(const FernDetection & fd) const
{
  FernDetection result = fd;
  result.featureData = new int[ivNumFerns];
  memcpy(result.featureData, fd.featureData, ivNumFerns*sizeof(int));
  return result;
}

#endif //FERNCLASSIFIER_H
