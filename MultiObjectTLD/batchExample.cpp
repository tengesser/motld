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

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#ifdef _MSC_VER
  #include "dirent.h"
#else
  #include <dirent.h>
#endif

#include "motld/MultiObjectTLD.h"
#include "motld/Utils.h"

#define DEFAULT_INPUT "input/motocross"
#define DEFAULT_OUTPUT "output"
#define MAX_FILE_NUMBER 0
#define OUTPUT_IMAGES 0

#define LOADCLASSIFIERATSTART 0
#define SAVECLASSIFIERATEND 0
#define LEARNMODEOFF 0
#define CLASSIFIERFILENAME "test.motld"

int pgm_select(const struct dirent *entry)
{
  return strstr(entry->d_name, ".pgm") != NULL;
}

int ppm_select(const struct dirent *entry)
{
  return strstr(entry->d_name, ".ppm") != NULL;
}

int main(int argc, char *argv[])
{
  std::string input_folder = DEFAULT_INPUT, output_folder = DEFAULT_OUTPUT;
  if(argc >= 2)
  {
    input_folder = argv[1];
    if(argc >= 3)
    {
      output_folder = argv[2];
    }
  }
  std::cout << "input folder: " << input_folder << std::endl;
  
  struct dirent **filelist;
  int fcount = -1;
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
    std::cout << "please create the file \"" << filename 
        << "\" specifying the initial bounding box[es] (x1,y1,x2,y2)" << std::endl;
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
  if (dir == NULL)
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
#if LOADCLASSIFIERATSTART
  MultiObjectTLD p = MultiObjectTLD::loadClassifier((char*)CLASSIFIERFILENAME);
#else
  MOTLDSettings settings(gray ? COLOR_MODE_GRAY : COLOR_MODE_RGB);
  MultiObjectTLD p(width, height, settings);
#endif
  
#if LEARNMODEOFF
  p.enableLearning(false);
#endif
  
  std::vector<ObjectBox> addBoxes;
  std::vector<ObjectBox>::iterator boxIt = boxes.begin();
  
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

    #if OUTPUT_IMAGES>0
    // and save debug image to file
    sprintf(filename, "%s/%s", output_folder.c_str(), filelist[i]->d_name);
    p.writeDebugImage(img,filename);
    #endif
    
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
  }
  outStream.close();

  std::cout << "MultiObjectTLD finished!" << std::endl;
#if SAVECLASSIFIERATEND
  std::cout << "Saving ..." << std::endl;
  p.saveClassifier((char*)CLASSIFIERFILENAME);
#endif
  for(int i = 0; i < fcount; i++)
    free(filelist[i]);
  free(filelist);
  
  return 0;
}
