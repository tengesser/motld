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

#ifndef LKTracker_H
#define LKTracker_H
  
#include "Matrix.h"
#include <vector>
#include <algorithm>
#include <math.h>

#define KERNEL_WIDTH 2
#define MAX_PYRAMID_LEVEL 5
/// number of iterations in each pyramid level
#define LK_ITERATIONS 30
/// number of tracking points in the uniform grid
#define GRID_SIZE_X 10
#define GRID_SIZE_Y 10
/// relative padding to borders of bounding box
/// @note original TLD uses an absolute padding of 5px
#define GRID_PADDING 0.15
/// number of (additional) points chosen by cornerness measure (0 = disable => more efficient)
#define N_CORNERNESS_POINTS 0
/// absolute padding (used for cornerness points only)
#define GRID_MARGIN 3
#define KERNEL_SIZE ((KERNEL_WIDTH*2+1)*(KERNEL_WIDTH*2+1))
  
/** @brief This class contains the "short term tracking" part of the algorithm.
 */
class LKTracker
{
public:
  /// Constructor
  LKTracker(int width, int height) : ivWidth(width), ivHeight(height), 
      ivPrevPyramid(NULL), ivIndex(1) {};
  /// Destructor
  ~LKTracker() {delete ivPrevPyramid;};
  /// Sets up the internal image pyramid
  void initFirstFrame(unsigned char * img);
  /// Sets up the internal image pyramid
  void initFirstFrame(const Matrix& img);
  /** @brief Computes the optical flow for each object
   *  @param bbox List of current object boxes, they are replaced with the new boxes
   *  @param isDefined Must have the same size as @b bbox. True for each object that is 
   *    currently defined and should be tracked. Is set to false if tracking failed.
   */
  void processFrame(const Matrix& curImage, std::vector<ObjectBox>& bbox, std::vector<bool>& isDefined); 
  /// An adapter for the single object case
  bool processFrame(const Matrix& curImage, ObjectBox& bbox, bool dotracking = true); 
  /// A list of points [x0,y0,...,xn,yn] that where considered as inliers in the last iteration
  const std::vector<int> * getDebugPoints() const { return &ivDebugPoints; };
  
private:
  /// Internal representation for an image pyramid
  struct LKPyramid
  {
    std::vector<Matrix> I,Ix,Iy;
    LKPyramid(){};
    LKPyramid(int nLevels){
      I = std::vector<Matrix>(nLevels);
      Ix = std::vector<Matrix>(nLevels);
      Iy = std::vector<Matrix>(nLevels); 
    };
  };
  /// Simple representation for 2D (sub pixel) image points
  struct Point2D
  {
      float x, y;
  };
  int ivWidth;
  int ivHeight;
  LKPyramid* ivPrevPyramid;
  int ivIndex;
  std::vector<int> ivDebugPoints;
  /** Computes median of a vector
   * @note changes order of vector-elements! */
  inline float median(std::vector<float> * vec, bool compSqrt = false) const;
  /** Computes normalized cross correlation
   * @details defined as: @f[NCC(A,B):=\frac{\sum_{x,y}(A(x,y)-\bar{A})(B(x,y)-\bar{B})}
    *    {\sqrt{\sum_{x,y}(A(x,y)-\bar{A})^2\sum_{x,y}(B(x,y)-\bar{B})^2}} @f]
   */
  inline double NCC(const Matrix& aMatrix, const Matrix& bMatrix) const;
  /** Computes optical flow for each tracking point. 
   * @details Based on the technical report "Pyramidal Implementation of the 
   *  Lucas Kanade Feature Tracker: Description of the algorithm" by Jean-Yves Bouguet */
  inline void pyramidLK(const LKPyramid *prevPyramid, const LKPyramid *curPyramid, 
                        const Point2D *prevPts, Point2D *nextPts, 
                        char *status, int count) const;
};



/**************************************************************************************************
 * IMPLEMENTATION                                                                                 *
 **************************************************************************************************/ 
void LKTracker::initFirstFrame(unsigned char * img)
{
  Matrix curImage(ivWidth, ivHeight);
  curImage.copyFromCharArray(img);
  initFirstFrame(curImage);
}
void LKTracker::initFirstFrame(const Matrix& img)
{
  ivPrevPyramid = new LKPyramid(MAX_PYRAMID_LEVEL+1);
  ivPrevPyramid->I[0] = img;
  for (int i = 0; i <= MAX_PYRAMID_LEVEL; ++i)
  {
    ivPrevPyramid->I[i].scharrDerivativeX(ivPrevPyramid->Ix[i]);
    ivPrevPyramid->I[i].scharrDerivativeY(ivPrevPyramid->Iy[i]);
    if (i < MAX_PYRAMID_LEVEL)
      ivPrevPyramid->I[i].halfSizeImage(ivPrevPyramid->I[i+1]);
    #if DEBUG > 1
    char filename[255];
    sprintf(filename, "output/img%05d-%d.ppm", 0, i);
    ivPrevPyramid->I[i].writeToPGM(filename);
    #endif    
  }
  #if DEBUG
  std::cout << "#1 LKTracker: initialized, image size = (" << img.xSize() << "," << img.ySize() << ")" << std::endl;
  #endif
}

bool LKTracker::processFrame(const Matrix& curImage, ObjectBox& bbox, bool dotracking)
{
  std::vector<ObjectBox> boxes;
  std::vector<bool> isDefined;
  boxes.push_back(bbox);
  isDefined.push_back(dotracking);
  processFrame(curImage, boxes, isDefined);
  bbox = boxes[0];
  return isDefined[0];
}

void LKTracker::processFrame(const Matrix& curImage, std::vector<ObjectBox>& bbox, std::vector<bool>& isDefined)
{
  int nobs = bbox.size();
  if (nobs > 0 && !ivPrevPyramid)
    initFirstFrame(curImage);
  #if DEBUG
  std::cout << "#" << (ivIndex+1) << " LKTracker: ";
  #endif
  ivDebugPoints.clear();
  LKPyramid* curPyramid = new LKPyramid(MAX_PYRAMID_LEVEL+1);
  curPyramid->I[0] = curImage;
  for (int i = 0; i < MAX_PYRAMID_LEVEL; ++i)
  {
    curPyramid->I[i].halfSizeImage(curPyramid->I[i+1]);
    #if DEBUG > 1
    char filename[255];
    sprintf(filename, "output/img%05d-%d.ppm", ivIndex, i);
    curPyramid->I[i].writeToPGM(filename);
    #endif  
  }
  
  #pragma omp parallel sections
  {
    #pragma omp section
    {
      //#pragma omp parallel for
      for (int i = 0; i <= MAX_PYRAMID_LEVEL; ++i)
        curPyramid->I[i].scharrDerivativeX(curPyramid->Ix[i]);
    }
    #pragma omp section
    {
      //#pragma omp parallel for
      for (int i = 0; i <= MAX_PYRAMID_LEVEL; ++i)
        curPyramid->I[i].scharrDerivativeY(curPyramid->Iy[i]);
    }
  }
  
  #if DEBUG > 1
  Matrix debugFlow(ivWidth, ivHeight, 0);
  #endif
      
  // loop over all object boxes
  for (int obj = 0; obj < nobs; obj++)
  {
    #if DEBUG
    std::cout << "\tObj" << obj << ": ";
    #endif
    if (isDefined[obj])
    {
      float oldwidth = bbox[obj].width, 
            oldheight = bbox[obj].height,
            oldcenterx = bbox[obj].x + oldwidth*0.5, 
            oldcentery = bbox[obj].y + oldheight*0.5;

      Point2D points0[N_CORNERNESS_POINTS + GRID_SIZE_X*GRID_SIZE_Y]; // points to be tracked
      Point2D points1[N_CORNERNESS_POINTS + GRID_SIZE_X*GRID_SIZE_Y]; // result of forward step
      Point2D points2[N_CORNERNESS_POINTS + GRID_SIZE_X*GRID_SIZE_Y]; // result of backward step
      char status[N_CORNERNESS_POINTS + GRID_SIZE_X*GRID_SIZE_Y];
      float fb[N_CORNERNESS_POINTS + GRID_SIZE_X*GRID_SIZE_Y];
      float ncc[N_CORNERNESS_POINTS + GRID_SIZE_X*GRID_SIZE_Y];
      int count = 0;
  
      // take points on a regular grid
      float stepX = oldwidth * (1 - 2*GRID_PADDING) / (GRID_SIZE_X-1), 
            stepY = oldheight * (1 - 2*GRID_PADDING) / (GRID_SIZE_Y-1);
      for (int y = 0; y < GRID_SIZE_Y; ++y)
        for (int x = 0; x < GRID_SIZE_X; ++x)
        {
          points0[count].x = bbox[obj].x + GRID_PADDING * oldwidth  + x*stepX;
          points0[count].y = bbox[obj].y + GRID_PADDING * oldheight + y*stepY;
          status[count] = 1;
          count++;
        }
    
      // ALTERNATIVE (or supplementary): compute interesting points (with high minimum eigen values) within bounding box
      #if N_CORNERNESS_POINTS > 0
      int xstart = std::max(0, (int)(bbox[obj].x)), 
          xend = std::min(ivWidth, (int)(bbox[obj].x + oldwidth + 1)),
          ystart = std::max(0, (int)(bbox[obj].y)),
          yend = std::min(ivHeight, (int)(bbox[obj].y + oldheight + 1)),
          xsize = xend - xstart,
          ysize = yend - ystart;      
      Matrix Ix2 (xsize, ysize);
      Matrix IxIy(xsize, ysize);
      Matrix Iy2 (xsize, ysize);
      Matrix cornerness(xsize, ysize, 0);
      for (int y = ystart; y < yend; y++)
        for (int x = xstart; x < xend; x++)
        {
          Ix2(x-xstart,y-ystart)  = ivPrevPyramid->Ix[0](x,y) * ivPrevPyramid->Ix[0](x,y);
          IxIy(x-xstart,y-ystart) = ivPrevPyramid->Ix[0](x,y) * ivPrevPyramid->Iy[0](x,y);
          Iy2(x-xstart,y-ystart)  = ivPrevPyramid->Iy[0](x,y) * ivPrevPyramid->Iy[0](x,y);
        }      
      Ix2.gaussianSmooth(2.0, 3);
      IxIy.gaussianSmooth(2.0, 3);
      Iy2.gaussianSmooth(2.0, 3);
      
      std::vector<float> cns;
      for (int y = GRID_MARGIN; y < ysize - GRID_MARGIN; y++)
        for (int x = GRID_MARGIN; x < xsize - GRID_MARGIN; x++)
        {
          cornerness(x,y) = (Ix2(x,y) + Iy2(x,y)) / 2.0 
                              - sqrt(((Ix2(x,y) + Iy2(x,y)) / 2.0) * ((Ix2(x,y) + Iy2(x,y)) / 2.0)
                                          - Ix2(x,y) * Iy2(x,y) + IxIy(x,y)*IxIy(x,y));
          if (cornerness(x,y) > 1.0)
            cns.push_back(cornerness(x, y));
        }
      float threshold = 0;
      if (cns.size() > N_CORNERNESS_POINTS)
      {
        std::nth_element(cns.begin(),cns.end() - N_CORNERNESS_POINTS, cns.end());
        threshold = *(cns.end() - N_CORNERNESS_POINTS);
      }
      for (int y = GRID_MARGIN; y < ysize - GRID_MARGIN; y++)
        for (int x = GRID_MARGIN; x < xsize - GRID_MARGIN; x++)
        {
          if (cornerness(x,y) > threshold && count < N_CORNERNESS_POINTS)
          {
            points0[count].x = x + xstart;
            points0[count].y = y + ystart;
            status[count] = 1;
            count++;
          }
        }      
      #endif //N_CORNERNESS_POINTS > 0
  
      // Track points forward
      pyramidLK(ivPrevPyramid, curPyramid, points0, points1, status, count);

      #if DEBUG > 1
      int nfwd = 0;
      for (int i = 0; i<count; ++i)
        if (status[i] > 0)
          ++nfwd;
      std::cout << "\t#fwd=" << nfwd;
      #endif
  
      // Track remaining points backward
      pyramidLK(curPyramid, ivPrevPyramid, points1, points2, status, count);
  
      // Compute FB-error and NCC
      std::vector<float> fbs, nccs;
      int nbwd = 0;
      //#pragma omp parallel for reduction(+:nbwd)
      for (int i = 0; i < count; ++i)
      {
        if (status[i] > 0)
        {
          fb[i] = sqrt((points2[i].x - points0[i].x) * (points2[i].x - points0[i].x) 
                     + (points2[i].y - points0[i].y) * (points2[i].y - points0[i].y));
          //#pragma omp critical
          fbs.push_back(fb[i]);
          Matrix mA = ivPrevPyramid->I[0].getRectSubPix(points0[i].x, points0[i].y, 10, 10);
          Matrix mB = curImage.getRectSubPix(points2[i].x, points2[i].y, 10, 10);
          ncc[i] = NCC(mA, mB);
          //#pragma omp critical
          nccs.push_back(ncc[i]);
          ++nbwd;
        }
      }
  
      float medFB = median(&fbs),
            medNCC = median(&nccs);
  
      #if DEBUG > 1
      std::cout << ", #bwd=" << nbwd;
      std::cout << "  \tmedFB=" << medFB << "\tmedNCC=" << medNCC;
      #endif
  
      //#pragma omp parallel for
      for (int i = 0; i < count; ++i)
      {
        if (status[i] > 0)
        {
          if (fb[i] > medFB || fb[i] > 8  || ncc[i] < medNCC)
            status[i] = 0;
          else 
          //#pragma omp critical
          {
            ivDebugPoints.push_back(round(points1[i].x)); 
            ivDebugPoints.push_back(round(points1[i].y)); 
          }
        }
      }
          
      // Compute median flow
      std::vector<float> deltax;
      std::vector<float> deltay;
      int num = 0;
      for (int i = 0; i < count; ++i)
      {
        if (status[i] > 0)
        {
          deltax.push_back(points1[i].x - points0[i].x);
          deltay.push_back(points1[i].y - points0[i].y);
          ++num;
          #if DEBUG > 1
          debugFlow.drawLine(points0[i].x, points0[i].y, points1[i].x, points1[i].y, 255);
          debugFlow.drawCross(points1[i].x, points1[i].y, 255);
          #endif
        }
      }
      if (num < 4)
      {
        #if DEBUG
        std::cout << "n=" << num << " => FAILURE: lost object" << std::endl;
        #endif
        isDefined[obj] = false;
        continue;
      }
      //else

      float dx = median(&deltax),
            dy = median(&deltay);  
      
      // Remove outliers
      /*
      for (int i = 0; i < count; ++i)
        if (status[i] > 0)
          if ((points1[i].x - points0[i].x - dx) * (points1[i].x - points0[i].x - dx)
              + (points1[i].y - points0[i].y - dy) * (points1[i].y - points0[i].y - dy) > 5*5)
          {
            status[i] = 0;
            num--;
          }
      */
  
      // Resize bounding box (compute median elongation factor)
      float s = 1;
      if (num >= 16){
        std::vector<float> d2;
        float dpx,dpy,ddx,ddy;
        for (int i = 0; i < count; ++i)
          if (status[i] > 0)
            for (int j = i + 1; j < count; ++j)
              if (status[j] > 0)
              {
                ddx = points0[i].x - points0[j].x;
                ddy = points0[i].y - points0[j].y;
                dpx = points1[i].x - points1[j].x;
                dpy = points1[i].y - points1[j].y;
                d2.push_back((dpx*dpx + dpy*dpy) / (ddx*ddx + ddy*ddy));
              }

        if (!d2.empty())
        {
          s = median(&d2, true);
          //upper bound for enlargement
          //s = std::min(1.1, s);
        }
      }
      //delete[] points0; delete[] points1; delete[] points2;
      //delete[] status; delete[] fb; delete[] ncc;

      float  centerx = oldcenterx + dx, 
            centery = oldcentery + dy;

      bbox[obj].x = (centerx - s * oldwidth * 0.5);
      bbox[obj].y = (centery - s * oldheight * 0.5);
      bbox[obj].width  = s * oldwidth;
      bbox[obj].height = s * oldheight;
      #if DEBUG
      std::cout << "n = " << num 
        << ", new BB: (" << round(bbox[obj].x) << "," << round(bbox[obj].y) << ", "
        << round(bbox[obj].width) << "," << round(bbox[obj].height) << ")";
        //<< ")\ts=" << s << std::endl;
      #endif  
    }else{
      #if DEBUG
      std::cout << "not defined";
      #endif
    }
  } // end for(obj)
  
  #if DEBUG > 1
  char filename[255];
  sprintf(filename, "output/flow%05d.ppm", ivIndex);
  writePPM(filename, ivPrevPyramid->I[0], curImage, debugFlow);
  #endif
  
  delete ivPrevPyramid;
  ivPrevPyramid = curPyramid;
  ++ivIndex;
}
 
inline void LKTracker::pyramidLK(const LKPyramid *prevPyramid, const LKPyramid *curPyramid, 
                          const Point2D *prevPts, Point2D *nextPts, 
                          char *status, int count) const
{
  for (int l = MAX_PYRAMID_LEVEL; l >= 0; --l)
  {
    int xSize = prevPyramid->I[l].xSize(), 
        ySize = prevPyramid->I[l].ySize();
    #if DEBUG > 2
    std::cout << "l=" << l << ", Size=(" << xSize << "," << ySize << ")" << std::endl;
    #endif
    
    #pragma omp parallel for default(shared)
    for (int i = 0; i < count; i++)
    {
      if (status[i] > 0)
      {
        //initial guess from previous iteration
        if (l == MAX_PYRAMID_LEVEL)
        {
          nextPts[i].x = prevPts[i].x;
          nextPts[i].y = prevPts[i].y;
        }else{
          nextPts[i].x *= 2.0;
          nextPts[i].y *= 2.0;
        }
        float px = prevPts[i].x * 1.0/(1<<l), py = prevPts[i].y * 1.0/(1<<l);
        int px0 = (int)px, py0 = (int)py;
        float pxa = px - px0, pya = py - py0;
        #if DEBUG > 2
        std::cout << "  p=(" << px0 << "+" << pxa << ", " << py0 << "+" << pya  << ")" << std::endl;
        #endif
        if (px < KERNEL_WIDTH || py < KERNEL_WIDTH || px >= xSize-KERNEL_WIDTH-1 
              || py >= ySize-KERNEL_WIDTH-1)
        {
          if (l >= MAX_PYRAMID_LEVEL-1){
            // Give it another try one level above
            nextPts[i].x = px;
            nextPts[i].y = py;
          }else{
            status[i] = 0;
          }
          //continue;        
        }else{ //omp parallel for does not like continues...
          // Compute components of spatial gradient Matrix G = [Gx2 Gxy; Gxy Gy2]
          float Gx2 = 0, Gxy = 0, Gy2 = 0;
          for (int x = px0-KERNEL_WIDTH; x <= px0+KERNEL_WIDTH+1; ++x)
          {
            float factor = (x == px0-KERNEL_WIDTH ? (1-pxa) : (x == px0+KERNEL_WIDTH+1 ? pxa : 1));
            for (int y = py0-KERNEL_WIDTH; y <= py0+KERNEL_WIDTH+1; ++y)
            {
              factor *= (y == py0-KERNEL_WIDTH ? (1-pya) : (y == py0+KERNEL_WIDTH+1 ? pya : 1));
              Gx2 += factor * prevPyramid->Ix[l](x,y)*prevPyramid->Ix[l](x,y); //Ix2(x,y)
              Gxy += factor * prevPyramid->Ix[l](x,y)*prevPyramid->Iy[l](x,y); //IxIy(x,y)
              Gy2 += factor * prevPyramid->Iy[l](x,y)*prevPyramid->Iy[l](x,y); //Iy2(x,y)
            }
          }    
          double denom = Gx2*Gy2 - Gxy*Gxy;
          #if DEBUG > 2
          std::cout << "\tGx2 = " << Gx2 << "\tGxy = " << Gxy << "\tGy2 = " << Gy2 << "\tdenom = " << denom << std::endl;
          #endif
          if (denom <= 1e-30)
          {
            status[i] = 0;
            //continue;        
          }else{          
            //iteratively compute additional flow on this pyramid level
            for (int k = 1; k <= LK_ITERATIONS; k++)
            {
              //float qx = px + tp->fx, qy = py + tp->fy;
              float qx = nextPts[i].x, qy = nextPts[i].y;
              if (qx < KERNEL_WIDTH || qy < KERNEL_WIDTH || qx >= xSize-KERNEL_WIDTH-1 || qy >= ySize-KERNEL_WIDTH-1)
              {  //lost tracking point
                if (l >= MAX_PYRAMID_LEVEL-1){
                  //give it another try one level above
                  nextPts[i].x = px;
                  nextPts[i].y = py;
                }else{
                  status[i] = 0;
                }
                break;
              }
              int vx0 = (int)qx - px0, vy0 = (int)qy - py0;
              float vxa = fmod(qx, 1), vya = fmod(qx, 1);
          
              //compute image missmatch vector b = [bx; by]
              float bx = 0, by = 0;
              for (int x = px0 - KERNEL_WIDTH; x <= px0 + KERNEL_WIDTH; ++x)
              {
                for (int y = py0 - KERNEL_WIDTH; y <= py0 + KERNEL_WIDTH; ++y)
                {
                  float dIk = (1-pxa) * ((1-pya)*prevPyramid->I[l](x, y) + pya*prevPyramid->I[l](x, y+1))
                              + pxa * ((1-pya)*prevPyramid->I[l](x+1, y) + pya*prevPyramid->I[l](x+1, y+1))
                              - (1-vxa) * ((1-vya)*curPyramid->I[l](x+vx0, y+vy0) + vya*curPyramid->I[l](x+vx0, y+vy0+1))
                              - vxa * ((1-vya)*curPyramid->I[l](x+vx0+1, y+vy0) + vya*curPyramid->I[l](x+vx0+1, y+vy0+1));
                  bx += dIk * ((1-pxa) * ((1-pya)*prevPyramid->Ix[l](x, y) + pya*prevPyramid->Ix[l](x, y+1)) 
                                  + pxa * ((1-pya)*prevPyramid->Ix[l](x+1, y) + pya*prevPyramid->Ix[l](x+1, y+1))); 
                  by += dIk * ((1-pxa) * ((1-pya)*prevPyramid->Iy[l](x, y) + pya*prevPyramid->Iy[l](x, y+1)) 
                                  + pxa * ((1-pya)*prevPyramid->Iy[l](x+1, y) + pya*prevPyramid->Iy[l](x+1, y+1))); 
                }
              }    
              float dx = (bx*Gy2 - by*Gxy) / denom,
                    dy = (by*Gx2 - bx*Gxy) / denom;
              nextPts[i].x += dx;
              nextPts[i].y += dy;
              #if DEBUG > 2
              std::cout << "\tf=(" << (nextPts[i].x - prevPts[i].x) << "," << (nextPts[i].y - prevPts[i].y) << ")" << std::endl;
              #endif
                        
              if (fabs(dx) > 3.5 || fabs(dy) > 3.5)
              {  //remove point because of unstable drifting..
                if (l >= 1){
                  nextPts[i].x = px;
                  nextPts[i].y = py;
                }else{
                  status[i] = 0;
                }
                break;
              }
            } //end for k
          }
        }
      } //end if (status > 0)
    } //end for each tp
  } //end for l
}

//compute median of a vector
//side effect: changes order of vector-elements!
inline float LKTracker::median(std::vector<float> * vec, bool compSqrt) const
{
  int n = vec->size();
  if (n == 0)
    return 0;
  if (n % 2) //odd: return (sqrt of) middle element
  { 
    std::nth_element(vec->begin(), vec->begin() + n/2, vec->end());
    return compSqrt ? sqrt(*(vec->begin() + n/2)) : *(vec->begin() + n/2);
  }else{     //even: return average (sqrt) of the two middle elements
    std::nth_element(vec->begin(), vec->begin() + (n/2-1), vec->end());
    float tmp = (compSqrt ? sqrt(*(vec->begin() + (n/2-1))) : *(vec->begin() + (n/2-1)));
    std::nth_element(vec->begin() + n/2, vec->begin() + n/2, vec->end());
    return 0.5 * (tmp + (compSqrt ? sqrt(*(vec->begin() + n/2)) : *(vec->begin() + n/2)));
  }
  return 0;
}

inline double LKTracker::NCC(const Matrix& aMatrix, const Matrix& bMatrix) const
{
  int size = aMatrix.size();
  if (size != bMatrix.size()){
    std::cerr << "ncc called for matrices with unequal size!" << std::endl;
    return 0;
  }    
  float aMean = aMatrix.avg(), bMean = bMatrix.avg();
  double sumA = 0, sumB = 0, sumDiff = 0;
  for (int i = 0; i < size; ++i)
  {
    sumA += (aMatrix.data()[i] - aMean) * (aMatrix.data()[i] - aMean);
    sumB += (bMatrix.data()[i] - bMean) * (bMatrix.data()[i] - bMean);
    sumDiff += (aMatrix.data()[i] - aMean) * (bMatrix.data()[i] - bMean);
  }
  if (sumA == 0 || sumB == 0)
    return 0;
  return sumDiff / sqrt(sumA*sumB); //+1)/2.0;
}

#endif 
