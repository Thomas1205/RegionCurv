/*** written by Thomas Schoenemann as an employee of Lund University, June 2010 ***/

#include "application.hh"
#include "grayimage.hh"
#include "colorimage.hh"
#include "lp_inpainting.hh"
#include "timing.hh"
#include "label_components.hh"
#include "color_conversion.hh"

int main(int argc, char** argv) {

  if (argc == 1 || (argc == 2 && strings_equal(argv[1],"-h"))) {

    std::cerr << "USAGE: " << argv[0] << std::endl
	      << "  -i <pgm> : filename of input image (to be inpainted)" << std::endl
	      << "  -mask <pgm> : filename of the image mask, values < 128 denote areas to be inpainted" << std::endl
	      << "  -lambda <double> : length weight" << std::endl
	      << "  -gamma <double> : curvature weight" << std::endl
	      << "  -o <filename> : name of the output inpainted image" << std::endl
	      << " [-n (4|8|16)]: size of neighborhood, default 8" << std::endl
	      << " [-b <uint>]: number of bins for the levels" << std::endl
	      << " [-enforce-levelcon]: enforce the level constraint" << std::endl
	      << "  -boundary-constraints (tight | simple | extra | off) : constraints to ensure consistency of regions and boundaries." 
	      << "     default is simple, extra unites simple and tight " << std::endl;
    exit(0);
  }

  const int nParams = 13;
  ParamDescr  params[nParams] = {{"-i",mandInFilename,0,""},{"-lambda",optWithValue,1,"0.001"},
				 {"-gamma",optWithValue,1,"1.0"},{"-o",mandOutFilename,0,""},
				 {"-mask",mandInFilename,0,""},{"-n",optWithValue,1,"8"},{"-b",optWithValue,1,"1"},
				 {"-boundary-constraints",optWithValue,1,"simple"},
				 {"-enforce-levelcon",flag,0,""},{"-light-constraints",flag,0,""},
				 {"-solver",optWithValue,1,"clp"},{"-curv-power",optWithValue,1,"2.0"},
				 {"-sequential",flag,0,""}};

  Application app(argc,argv,params,nParams);

  Math3D::NamedColorImage<float> color_image(app.getParam("-i"),MAKENAME(color_image));
  Math3D::NamedColorImage<float> yuv_image(MAKENAME(yuv_image));

  uint xDim = color_image.xDim();
  uint yDim = color_image.yDim();
  uint zDim = color_image.zDim();

  Math3D::NamedColorImage<float> yuv_out_image(xDim,yDim,zDim,MAKENAME(yuv_out_image));
  Math3D::NamedColorImage<float> rgb_out_image(xDim,yDim,zDim,MAKENAME(rgb_out_image));

  bool enforce_boundary_consistency = false;
  bool enforce_region_edge_consistency = false;

  std::string constraint_string = app.getParam("-boundary-constraints");
  if (constraint_string == "tight") {
    enforce_region_edge_consistency = true;
  }
  else if (constraint_string == "simple") {
    enforce_boundary_consistency = true;
  }
  else if (constraint_string == "extra") {
    enforce_region_edge_consistency = true;
    enforce_boundary_consistency = true;
  }
  else if (constraint_string != "off") {
    USER_ERROR << "unknown boundary constraint mode \"" << constraint_string << "\"" << std::endl
	       << " choose from (tight | simple | extra | off)" << std::endl
	       << "  Exiting." << std::endl;
    exit(1);
  }


  if (zDim == 1)
    yuv_image = color_image;
  else {
    yuv_image.resize(xDim,yDim,3);

    uchar yy,uu,vv;
    for (uint y=0; y < yDim; y++) {
      for (uint x=0; x < xDim; x++) {
	rgb2yuv(color_image(x,y,0),color_image(x,y,1),color_image(x,y,2),yy,uu,vv);
	yuv_image(x,y,0) = yy;
	yuv_image(x,y,1) = uu;
	yuv_image(x,y,2) = vv;
      }
    }
  }

  Math2D::NamedGrayImage<float> mask(app.getParam("-mask"),MAKENAME(mask));

  if (mask.xDim() != xDim || mask.yDim() != yDim) {
    USER_ERROR << "dimensions of image and mask mismatch" << std::endl;
    exit(1);
  }
  
  Math2D::Matrix<uint> component;
  uint nComponents = label_components(mask,component);

  double lambda = convert<double>(app.getParam("-lambda"));
  double gamma = convert<double>(app.getParam("-gamma"));
  double curv_power = convert<double>(app.getParam("-curv-power"));

  std::clock_t tStartComputation, tEndComputation;
  tStartComputation = std::clock();

  Math2D::NamedMatrix<float> inpainted_image(xDim,yDim,0.0,MAKENAME(inpainted_image));

  uint neighborhood = convert<uint>(app.getParam("-n"));
  uint nBins = convert<uint>(app.getParam("-b"));

  for (uint z=0; z < yuv_image.zDim(); z++) {

    Math2D::NamedMatrix<float> image(xDim,yDim,0.0,MAKENAME(image));
    for (uint y=0; y < yDim; y++) 
      for (uint x=0; x < xDim; x++)
	image(x,y) = yuv_image(x,y,z);


    float min_val = 1e36;
    float max_val = 0.0;
    
    
    for (uint y=0; y < yDim; y++) {
      for (uint x=0; x < xDim; x++) {
	if (mask(x,y) >= 128.0) {
	  
	  if (min_val > image(x,y))
	    min_val = image(x,y);
	  if (max_val < image(x,y))
	    max_val = image(x,y);
	}
      }
    }
    
    inpainted_image = image;

    for (uint y=0; y < yDim; y++) {
      for (uint x=0; x < xDim; x++) {
	if (mask(x,y) < 128.0) {
	  inpainted_image(x,y) = min_val;
	}
      }
    }
    

    for (uint c=1; c <= nComponents; c++) {
      
      std::cerr << "###################### channel " << (z+1) << "/" << zDim << ", component " << c << "/" << nComponents << std::endl;
      
      for (uint y=0; y < yDim; y++)
	for (uint x=0; x < xDim; x++)
	  mask(x,y) = (component(x,y) == c) ? 0.0 : 255.0;

      if (app.is_set("-sequential")) {
	
	bool first_level_found = false;
	
	for (uint level = uint(min_val)+1; level <= uint(max_val); level++) {
	  
	  std::cerr << "processing level: " << level << std::endl;
	  
	  Math2D::NamedMatrix<float> temp_image(xDim,yDim,0.0,MAKENAME(temp_image));
	  
	  for (uint y=0; y < yDim; y++) {
	    for (uint x=0; x < xDim; x++) {
	      if (image(x,y) >= level)
		temp_image(x,y) = 0;
	      else
		temp_image(x,y) = 1;
	    }
	  }
	  
	  double energy = lp_inpaint(temp_image, mask, lambda, gamma, curv_power, neighborhood, 0.0, app.getParam("-solver"),
				     temp_image, enforce_boundary_consistency, enforce_region_edge_consistency,
				     app.is_set("-light-constraints"));
	  
	  if (energy != 0.0)
	    first_level_found = true;
	  if (energy == 0.0 && first_level_found)
	    break;
	  
	  for (uint y=0; y < yDim; y++) {
	    for (uint x=0; x < xDim; x++) {
	      if (component(x,y) == c)
		inpainted_image(x,y) += 1.0 - temp_image(x,y);
	    }
	  }
	}
      }
      else if (nBins == 1) {
	lp_inpaint(inpainted_image, mask, lambda, gamma, curv_power, neighborhood, 0.0, app.getParam("-solver"),
		   inpainted_image, enforce_boundary_consistency, enforce_region_edge_consistency,
		   app.is_set("-light-constraints"));
      }
      else {
	
	lp_inpaint_hybrid(inpainted_image, mask, lambda, gamma, curv_power, neighborhood,nBins, 0.0, app.getParam("-solver"),
			  inpainted_image, enforce_boundary_consistency, enforce_region_edge_consistency,
			  app.is_set("-enforce-levelcon"),app.is_set("-light-constraints"));
      }
    }

    
    for (uint y=0; y < yDim; y++) 
      for (uint x=0; x < xDim; x++)
	yuv_out_image(x,y,z) = inpainted_image(x,y);
  }
    
  if (zDim == 3) {
    //back-transform
    uchar r,g,b;

    for (uint y=0; y < yDim; y++) {
      for (uint x=0; x < xDim; x++) {
	yuv2rgba(yuv_image(x,y,0),yuv_out_image(x,y,1),yuv_out_image(x,y,2),r,g,b);

	rgb_out_image(x,y,0) = r;
	rgb_out_image(x,y,1) = g;
	rgb_out_image(x,y,2) = b;
      }
    }
  }
  else
    rgb_out_image = yuv_out_image;

  tEndComputation = std::clock();

  rgb_out_image.savePPM(app.getParam("-o"),rgb_out_image.max_intensity());

  std::cerr << "computation took " << diff_seconds(tEndComputation, tStartComputation) << " seconds." << std::endl;

}
