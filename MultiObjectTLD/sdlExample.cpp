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
To run this example you need the libraries SDL and SDL_image

In Linux you can install them via packet manager (libsdl1.2-dev and libsdl-image1.2-dev)
or download and compile both from source (.configure / make / sudo make install)

For Windows you have to download the development versions and 
link the include-folders, the .lib's, and copy the .dll's
*/

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS 1
#pragma warning(disable: 4244) // possible loss of data
#endif

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#ifdef _MSC_VER
  #include "dirent.h"
#else
  #include <dirent.h>
#endif
#include <SDL.h>
#include <SDL_image.h>
#include <SDL_thread.h>

#include "motld/MultiObjectTLD.h"
#include "motld/Utils.h"

#define DEFAULT_INPUT "input/motocross"
#define DEFAULT_OUTPUT "output"
#define MAX_FILE_NUMBER 0

int pgm_select(const struct dirent *entry)
{
  return strstr(entry->d_name, ".pgm") != NULL;
}

int ppm_select(const struct dirent *entry)
{
  return strstr(entry->d_name, ".ppm") != NULL;
}

class SDLEngine
{
public:
  SDLEngine(){};
  ~SDLEngine(){SDL_Quit();};
  static void Init(int argc, char *argv[]);
  static int Run(void*);
  static bool ivQuit;
  static void Render();
  static void HandleInput();
private:
  static void SetSize(const int w, const int h);
  static SDL_Surface* ivScreen;
  static std::string input_folder, output_folder;
  static struct dirent **filelist;
  static int fcount;
  static int display, displaymax, curdisplay;
  static bool displaylast;
};

int main(int argc, char *argv[])
{
  SDLEngine::Init(argc, argv);
  
  SDL_CreateThread(SDLEngine::Run, 0);
  //engine.Run();
  while(!SDLEngine::ivQuit){
    SDLEngine::HandleInput();
    SDLEngine::Render();
    SDL_Delay(50);
  }
  SDL_Quit();
  return 0;
}

SDL_Surface* SDLEngine::ivScreen;
std::string SDLEngine::input_folder, SDLEngine::output_folder;
struct dirent **SDLEngine::filelist;
int SDLEngine::fcount = -1;
int SDLEngine::display = 0, SDLEngine::curdisplay = -1, SDLEngine::displaymax = -1;
bool SDLEngine::displaylast = true, SDLEngine::ivQuit = false;

void SDLEngine::Init(int argc, char *argv[])
{
  atexit( SDL_Quit );
  if(SDL_Init(SDL_INIT_VIDEO) < 0)
    std::cout << "SDLerror: " << SDL_GetError() << std::endl;
  input_folder = DEFAULT_INPUT; output_folder = DEFAULT_OUTPUT;
  if(argc >= 2)
  {
    input_folder = argv[1];
    if(argc >= 3)
    {
      output_folder = argv[2];
    }
  }
  SDL_EnableKeyRepeat(30, 30);
  ivScreen = SDL_SetVideoMode( 400, 300, 0, SDL_HWSURFACE | SDL_DOUBLEBUF); // | SDL_RESIZABLE
  if(ivScreen == NULL)
    std::cout << "ERROR: " << SDL_GetError() << std::endl;
  SDL_WM_SetCaption("MultiObjectTLD", 0 );
}

int SDLEngine::Run(void*)
{
  std::cout << "input folder: " << input_folder << std::endl;
  bool gray = false;
  fcount = scandir(input_folder.c_str(), &filelist, ppm_select, alphasort);
  if (fcount <= 0)
  {
    fcount = scandir(input_folder.c_str(), &filelist, pgm_select, alphasort);
    gray = true;
  }
  if (fcount <= 0)
  {
    std::cout << "There are no .ppm or .pgm files in this folder! Maybe you have to convert the images first e.g. using" << std::endl
      << "  mogrify -format ppm *.jpg" << std::endl;
    return 0;
  }
  std::cout << "found " << fcount << " files" << std::endl;
  char filename[255];
  sprintf(filename, "%s/init.txt", input_folder.c_str());
  std::ifstream aStream(filename);
  if(!aStream || aStream.eof())
  {
    std::cout << "please create the file \"" << filename << "\" specifying the initial bounding box (x1,y1,x2,y2)" << std::endl;
    return 0;
  }
  char line[255];
  int x1,y1,x2,y2,imgid, width,height;
  std::vector<ObjectBox> boxes;
  while(aStream.getline(line,255))
  {
    x1 = y1 = x2 = y2 = imgid = 0;
    int i = 0;
    for(;line[i] >= '0' && line[i] <= '9'; i++)
      x1 = x1*10 + (line[i] - '0');
    for(i++;line[i] >= '0' && line[i] <= '9'; i++)
      y1 = y1*10 + (line[i] - '0');
    for(i++;line[i] >= '0' && line[i] <= '9'; i++)
      x2 = x2*10 + (line[i] - '0');
    for(i++;line[i] >= '0' && line[i] <= '9'; i++)
      y2 = y2*10 + (line[i] - '0');
    if(line[i] == ',')
      for(i++;line[i] >= '0' && line[i] <= '9'; i++)
        imgid = imgid*10 + (line[i] - '0'); 
    ObjectBox b = {(float)x1,(float)y1,(float)(x2-x1),(float)(y2-y1),imgid};
    boxes.push_back(b);
  }
  aStream.close();
    
  std::cout << "output folder: " << output_folder << std::endl;
  DIR * dir = opendir(output_folder.c_str());
  if (dir == 0)
  {
    std::cout << "\tdoes not exist -> try to create it" << std::endl;
    if(system(("mkdir "+output_folder).c_str()))
    {
      std::cout << "\t failed to create directory" << std::endl;
      return 0;
    }
  }
  closedir(dir);

  sprintf(filename, "%s/%s", input_folder.c_str(), filelist[0]->d_name);
  int z;
  unsigned char* dummy = gray ? readFromPGM<unsigned char>(filename, width, height) :
                                readFromPPM<unsigned char>(filename, width, height, z);
  delete[] dummy;
    
  // Initialize MultiObjectTLD
  MOTLDSettings s(gray ? COLOR_MODE_GRAY : COLOR_MODE_RGB);
  // s.bbMin = 18;
  MultiObjectTLD p(width, height, s);
  std::vector<ObjectBox> addBoxes;
  std::vector<ObjectBox>::iterator boxIt = boxes.begin();

  if(ivScreen != NULL)
    SDL_FreeSurface(ivScreen);
  ivScreen = SDL_SetVideoMode( width, height, 0, SDL_HWSURFACE | SDL_DOUBLEBUF); // | SDL_RESIZABLE
  SDL_WM_SetCaption("MultiObjectTLD", 0 );

  sprintf(filename, "%s/output.txt", output_folder.c_str());
  std::ofstream outStream(filename);  

  for (int i=0; i < fcount && (!MAX_FILE_NUMBER || i<MAX_FILE_NUMBER); ++i)
  {
    // first load the image
    sprintf(filename, "%s/%s", input_folder.c_str(), filelist[i]->d_name);
    int xS, yS, z;
    unsigned char* img = gray ? readFromPGM<unsigned char>(filename, xS, yS) :
                                readFromPPM<unsigned char>(filename, xS, yS, z);
    // then process it with MultiObjectTLD
    p.processFrame(img);
    
    while(boxIt != boxes.end() && boxIt->objectId == i)
    {
      addBoxes.push_back(*boxIt);
      boxIt++;
    }
    if(addBoxes.size() > 0){
      p.addObjects(addBoxes);
      addBoxes.clear();
    }
    
    // and save debug image to file
    sprintf(filename, "%s/%s", output_folder.c_str(), filelist[i]->d_name);
    p.writeDebugImage(img,filename);
    displaymax = i;

    // print current box to output file
    if(p.getValid())
    {
      ObjectBox b = p.getObjectBox();
      if(i > 0)
        outStream << std::endl;
      outStream << b.x << "," << b.y << "," << (b.x + b.width) << "," << (b.y + b.height);
    }
    else
      outStream << std::endl << "NaN,NaN,NaN,NaN";
    delete[] img;
    if(ivQuit)
      break;
  }
  outStream.close();
  std::cout << "MultiObjectTLD finished!" << std::endl;
  return 0;
}

void SDLEngine::Render()
{
  if(displaymax < 0)
    return;
  //SDL_FillRect( ivScreen, 0, SDL_MapRGB( ivScreen->format, 255, 255, 255 ) );
  if(displaylast)
    display = displaymax;
  char filename[255];
  if(display != curdisplay)
  {
    sprintf(filename, "%s/%s", output_folder.c_str(), filelist[display]->d_name);
    SDL_Surface* sdlimg = IMG_Load(filename);
    SDL_Rect rect = {0,0,(Uint16)ivScreen->w,(Uint16)ivScreen->h};
    SDL_BlitSurface(sdlimg, &rect, ivScreen, &rect);
    SDL_Flip(ivScreen);
    SDL_FreeSurface(sdlimg);
    curdisplay = display;
  }
  sprintf(filename, "MultiObjectTLD: %d / %d", display+1, displaymax+1);
  SDL_WM_SetCaption(filename, 0);
}

void SDLEngine::HandleInput()
{
  // Poll for events, and handle the ones we care about.
  SDL_Event event;
  while ( SDL_PollEvent( &event ) ) 
  {
    switch ( event.type ) 
    {
    case SDL_KEYDOWN:
      switch (event.key.keysym.sym)
      {
      case SDLK_LEFT:
        displaylast = false;
        display--;
        if(display < 0)
          display = 0;
        break;
      case SDLK_RIGHT:
        displaylast = false;
        display++;
        if(display > displaymax)
          display = displaymax;        
        break;
      case SDLK_PAGEUP:
        displaylast = false;
        display -= 10;
        if(display < 0)
          display = 0;
        break;
      case SDLK_PAGEDOWN:
        displaylast = false;
        display += 10;
        if(display > displaymax)
          display = displaymax;        
        break;
      case SDLK_END:
        displaylast = true;
        display = displaymax;
        break;
      case SDLK_HOME:
        displaylast = false;
        display = 0;
        break;
      case SDLK_ESCAPE:  
        ivQuit = true;
        break;
      default: 
        break;
      }
       break;

    case SDL_QUIT:
      ivQuit = true;
      break;

    } // switch
  } // while (handling input)
}
 
