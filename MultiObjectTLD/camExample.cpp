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
 
/*
This is a live demo invoking a webcam.
It makes use of OpenCV to capture the frames and highgui to display the output.

Use your mouse to draw bounding boxes for tracking.
There are some keys to customize which components are displayed:
 D - dis/enable drawing of detections (green boxes)
 P - dis/enable drawing of learned patches
 T - dis/enable drawing of tracked points
 L - dis/enable learning
 S - save current classifier to CLASSIFIERFILENAME
 O - load classifier from CLASSIFIERFILENAME (not implemented)
 R - reset classifier (not implemented)
 ESC - exit
*/

#include <stdlib.h>
#include <stdio.h>

#include "cv.h" 
#include "highgui.h" 

#include "motld/MultiObjectTLD.h"

#define LOADCLASSIFIERATSTART 0
#define CLASSIFIERFILENAME "test.moctld"

//uncomment if you have a high resolution camera and want to speed up tracking
#define FORCE_RESIZING
#define RESOLUTION_X 320
#define RESOLUTION_Y 240

#define MOUSE_MODE_MARKER 0
#define MOUSE_MODE_ADD_BOX 1
#define MOUSE_MODE_IDLE 2
IplImage* curImage = NULL;
bool ivQuit = false;
int ivWidth, ivHeight;
CvCapture* capture;
ObjectBox mouseBox = {0,0,0,0,0};
int mouseMode = MOUSE_MODE_IDLE;
int drawMode = 255;
bool learningEnabled = true, save = false, load = false, reset = false;

void Init(int argc, char *argv[]);
void* Run(void*);
void HandleInput(int interval = 1);
void MouseHandler(int event, int x, int y, int flags, void* param);
void FromRGB(Matrix& maRed, Matrix& maGreen, Matrix& maBlue);

int main(int argc, char *argv[])
{
  Init(argc, argv);
  Run(0);
  cvDestroyAllWindows();
  return 0;
}


void Init(int argc, char *argv[])
{  
  capture = cvCaptureFromCAM(CV_CAP_ANY);
  if(!capture){
    std::cout << "error starting video capture" << std::endl;
    exit(0);
  }
  //propose a resolution
  cvSetCaptureProperty(capture, CV_CAP_PROP_FRAME_WIDTH, 640);
  cvSetCaptureProperty(capture, CV_CAP_PROP_FRAME_HEIGHT, 480);
  //get the actual (supported) resolution
  ivWidth = cvGetCaptureProperty(capture, CV_CAP_PROP_FRAME_WIDTH);
  ivHeight = cvGetCaptureProperty(capture, CV_CAP_PROP_FRAME_HEIGHT);
  std::cout << "camera/video resolution: " << ivWidth << "x" << ivHeight << std::endl;
  #ifdef FORCE_RESIZING
  ivWidth = RESOLUTION_X;
  ivHeight = RESOLUTION_Y;
  #endif
  
  cvNamedWindow("MOCTLD", 0); //CV_WINDOW_AUTOSIZE );
  
  CvSize wsize = {ivWidth, ivHeight};
  curImage = cvCreateImage(wsize, IPL_DEPTH_8U, 3);
  
  cvResizeWindow("MOCTLD", ivWidth, ivHeight);
  cvSetMouseCallback("MOCTLD", MouseHandler);
}

void* Run(void*)
{
  int size = ivWidth*ivHeight;
    
  // Initialize MultiObjectTLD
  #if LOADCLASSIFIERATSTART
  MultiObjectTLD p = MultiObjectTLD::loadClassifier((char*)CLASSIFIERFILENAME);
  #else
  MOTLDSettings settings(COLOR_MODE_RGB);
  settings.useColor = true;
  MultiObjectTLD p(ivWidth, ivHeight, settings);
  #endif
  
  Matrix maRed;
  Matrix maGreen;
  Matrix maBlue;
  unsigned char img[size*3];
  #ifdef FORCE_RESIZING
  CvSize wsize = {ivWidth, ivHeight};
  IplImage* frame = cvCreateImage(wsize, IPL_DEPTH_8U, 3);
  #endif
  while(!ivQuit)
  {    
    /*
    if(reset){
      p = *(new MultiObjectTLD(ivWidth, ivHeight, COLOR_MODE_RGB));
      reset = false;
    }
    if(load){
      p = MultiObjectTLD::loadClassifier(CLASSIFIERFILENAME);
      load = false;
    }
    */
    
    // Grab an image
    if(!cvGrabFrame(capture)){
      std::cout << "error grabbing frame" << std::endl;
      break;
    }
    #ifdef FORCE_RESIZING
    IplImage* capframe = cvRetrieveFrame(capture);
    cvResize(capframe, frame);
    #else
    IplImage* frame = cvRetrieveFrame(capture);
    #endif
    for(int j = 0; j<size; j++){
      img[j] = frame->imageData[j*3+2];
      img[j+size] = frame->imageData[j*3+1];
      img[j+2*size] = frame->imageData[j*3];
    }
    
    // Process it with motld
    p.processFrame(img);
    
    // Add new box
    if(mouseMode == MOUSE_MODE_ADD_BOX){
      p.addObject(mouseBox);
      mouseMode = MOUSE_MODE_IDLE;
    }
    
    // Display result
    HandleInput();
    p.getDebugImage(img, maRed, maGreen, maBlue, drawMode);    
    FromRGB(maRed, maGreen, maBlue);
    cvShowImage("MOCTLD", curImage);
    p.enableLearning(learningEnabled);
    if(save){
      p.saveClassifier((char*)CLASSIFIERFILENAME);
      save = false;
    }
  }
  //delete[] img;
  cvReleaseCapture(&capture);
  return 0;
}

void HandleInput(int interval)
{
  int key = cvWaitKey(interval);
  if(key >= 0)
  {
    switch (key)
    {
      case 'd': drawMode ^= DEBUG_DRAW_DETECTIONS;  break;
      case 't': drawMode ^= DEBUG_DRAW_CROSSES;  break;
      case 'p': drawMode ^= DEBUG_DRAW_PATCHES;  break;
      case 'l':
        learningEnabled = !learningEnabled;
        std::cout << "learning " << (learningEnabled? "en" : "dis") << "abled" << std::endl;
        break;
      case 'r': reset = true; break;
      case 's': save = true;  break;
      case 'o': load = true;  break;
      case 27:  ivQuit = true; break; //ESC
      default: 
        //std::cout << "unhandled key-code: " << key << std::endl;
        break;
    }
  }
}

void MouseHandler(int event, int x, int y, int flags, void* param)
{
  switch(event){
    case CV_EVENT_LBUTTONDOWN:
      mouseBox.x = x;
      mouseBox.y = y;
      mouseBox.width = mouseBox.height = 0;
      mouseMode = MOUSE_MODE_MARKER;
      break;
    case CV_EVENT_MOUSEMOVE:
      if(mouseMode == MOUSE_MODE_MARKER){
        mouseBox.width = x - mouseBox.x;
        mouseBox.height = y - mouseBox.y;
      }
      break;
    case CV_EVENT_LBUTTONUP:
      if(mouseMode != MOUSE_MODE_MARKER)
        break;
      if(mouseBox.width < 0){
        mouseBox.x += mouseBox.width;
        mouseBox.width *= -1;
      }
      if(mouseBox.height < 0){
        mouseBox.y += mouseBox.height;
        mouseBox.height *= -1;
      }
      if(mouseBox.width < 4 || mouseBox.height < 4){
        std::cout << "bounding box too small!" << std::endl;
        mouseMode = MOUSE_MODE_IDLE;
      }else
        mouseMode = MOUSE_MODE_ADD_BOX;
      break;
    case CV_EVENT_RBUTTONDOWN:
      mouseMode = MOUSE_MODE_IDLE;
      break;
  }
}

void FromRGB(Matrix& maRed, Matrix& maGreen, Matrix& maBlue)
{
  for(int i = 0; i < ivWidth*ivHeight; ++i){
    curImage->imageData[3*i+2] = maRed.data()[i];
    curImage->imageData[3*i+1] = maGreen.data()[i];
    curImage->imageData[3*i+0] = maBlue.data()[i];
  }
  //at this place you could save the images using
  //cvSaveImage(filename, curImage);
  if(mouseMode == MOUSE_MODE_MARKER)
  {
    CvPoint pt1; pt1.x = mouseBox.x; pt1.y = mouseBox.y;
    CvPoint pt2; pt2.x = mouseBox.x + mouseBox.width; pt2.y = mouseBox.y + mouseBox.height;  
    cvRectangle(curImage, pt1, pt2, CV_RGB(0,0,255));
  }
}
