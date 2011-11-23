/*** written by Thomas Schoenemann as an employee of Lund University, Sweden, September 2010 ***/

#include "application.hh"
#include "grayimage.hh"
#include "colorimage.hh"
#include "conversion.hh"
#include "draw_segmentation.hh"
#include "lp_segmentation.hh"
#include "timing.hh"
#include "smoothing.hh"

void check_filename(std::string name)
{
  std::ofstream fout(name.c_str());
  if (!fout) {
    std::cerr << "Cannot open " << name << std::endl;
    exit(2);
  }
}

int main(int argc, char** argv) {

  if (argc == 1 || (argc == 2 && strings_equal(argv[1],"-h"))) {

    std::cerr << "USAGE: " << argv[0] << std::endl
	      << "  -i <ppm or ppm> : filename of input image (to be segmented)" << std::endl
	      << "  -lambda <double> : length weight" << std::endl
	      << "  -gamma <double> : curvature weight" << std::endl
	      << "  -o <filename> : name of the output segmentation" << std::endl
	      << "  -fg-mask <filename> : file containing foreground seeds" << std::endl
	      << "  -bg-mask <filename> : file containing background seeds" << std::endl
	      << "  -boundary-constraints (tight | simple | extra | off) : constraints to ensure consistency of regions and boundaries." 
	      << "     default is tight (= Strandmark&Kahl 11), extra unites simple and tight " << std::endl 
	      << " [-n (4|8|16)]: size of neighborhood, default 8" << std::endl
	      << " [-bruckstein]: curvature discretization according to Bruckstein et al." << std::endl
	      << " [-ignore-crossings] : allow crossings of line pairs, e.g. when three regions meet in a point" << std::endl
              << " [-light-constraints] : take only half of the optional constraints" << std::endl
	      << " [-reduce-pairs] : same some memory and run-time by not considering pairs with very high curvature" << std::endl
	      << " -solver ( clp | gurobi | mosek | cplex | xpress | own-conv ) : default clp" << std::endl;

    exit(0);
  }

 
  ParamDescr  params[] = {{"-i",mandInFilename,0,""},{"-lambda",optWithValue,1,"1.0"},
                          {"-gamma",optWithValue,1,"1.0"},{"-o",mandOutFilename,0,""},
                          {"-n",optWithValue,1,"8"},{"-boundary-constraints",optWithValue,0,"tight"},{"-light-constraints",flag,0,""},
                          {"-bruckstein",flag,0,""},{"-solver",optWithValue,1,"clp"},
                          {"-fg-mask",mandInFilename,0,""},{"-bg-mask",mandInFilename,0,""},
                          {"-griddim",optWithValue,1,"-1"},{"-debug-svg",flag,0,""},{"-ignore-crossings",flag,0,""},
			  {"-no-touching-regions",flag,0,""},{"-reduce-pairs",flag,0,""}};

  const int nParams = sizeof(params)/sizeof(ParamDescr);

  Application app(argc,argv,params,nParams);

  Math3D::NamedColorImage<float> color_image(app.getParam("-i"),MAKENAME(color_image));

  std::string base_filename = app.getParam("-o");
  //check_filename(base_filename + ".final.svg");
  //check_filename(base_filename + ".lp.svg");
  //check_filename(base_filename + ".lp_simple.svg");
  check_filename(base_filename + ".out.ppm");
  check_filename(base_filename + ".seg.pgm");
  //check_filename(base_filename + ".frac.pgm");

  uint xDim = uint( color_image.xDim() );
  uint yDim = uint( color_image.yDim() );
  uint zDim = uint( color_image.zDim() );

  Math2D::NamedGrayImage<float> fg_mask(app.getParam("-fg-mask"),MAKENAME(fg_mask));
  Math2D::NamedGrayImage<float> bg_mask(app.getParam("-bg-mask"),MAKENAME(bg_mask));

  if (fg_mask.xDim() != xDim || fg_mask.yDim() != yDim) {
    USER_ERROR << " dimensions of foreground mask do not match the supplied image" << std::endl;
    exit(1);
  }
  
  if (bg_mask.xDim() != xDim || bg_mask.yDim() != yDim) {
    USER_ERROR << " dimensions of background mask do not match the supplied image" << std::endl;
    exit(1);
  }

  uint nBins = 16;
  uint hist_dim = 256 / nBins;

  Math3D::Tensor<float> fg_hist(hist_dim,hist_dim,hist_dim,1e-2f);
  Math3D::Tensor<float> bg_hist(hist_dim,hist_dim,hist_dim,1e-2f);

  /*** set seeds and compute histograms ***/
  Math2D::NamedMatrix<int> seeds(xDim,yDim,-1,MAKENAME(seeds));

  for (uint y=0; y < yDim; y++) {
    for (uint x=0; x < xDim; x++) {

      uint r = uint( color_image(x,y,0) );
      uint g = 0;
      if (zDim > 1)
	g = uint( color_image(x,y,1) );
      uint b = 0;
      if (zDim > 1)
	b = uint( color_image(x,y,2) );

      if (fg_mask(x,y) < 127.0) {
	seeds(x,y) = 1;
	fg_hist(r/nBins, b/nBins, g/nBins) += 1.0;
      }
      else if (bg_mask(x,y) < 127.0) {
	seeds(x,y) = 0;
	bg_hist(r/nBins, b/nBins, g/nBins) += 1.0;
      }
    }
  }

  fg_hist *= 1.0f / fg_hist.sum();
  bg_hist *= 1.0f / bg_hist.sum();

  //smooth the histograms
  vol_smooth_isotropic_gauss(fg_hist, 0.75);
  vol_smooth_isotropic_gauss(bg_hist, 0.75);

  fg_hist *= 1.0 / fg_hist.sum();
  bg_hist *= 1.0 / bg_hist.sum();

  double lambda = convert<double>(app.getParam("-lambda"));
  double gamma = convert<double>(app.getParam("-gamma"));

  Math2D::NamedMatrix<float> data_term(xDim,yDim,MAKENAME(data_term));
  Math2D::NamedMatrix<uint> segmentation(xDim,yDim,0,MAKENAME(segmentation));

  double energy_offset = 0.0;

  for (uint y=0; y < yDim; y++) {
    for (uint x=0; x < xDim; x++) {

      uint r = uint( color_image(x,y,0) );
      uint g = 0;
      if (zDim > 1)
	g = uint(color_image(x,y,1));
      uint b = 0;
      if (zDim > 1)
	b = uint(color_image(x,y,2));

      double fg_term = - std::log(fg_hist(r/nBins, b/nBins, g/nBins)); 
      double bg_term = - std::log(bg_hist(r/nBins, b/nBins, g/nBins));

      data_term(x,y) = float(fg_term - bg_term);

      energy_offset += bg_term;
    }
  }

  std::clock_t tStartComputation, tEndComputation;
  tStartComputation = std::clock();

  LPSegOptions seg_opts;
  seg_opts.neighborhood_ = convert<uint>(app.getParam("-n"));
  seg_opts.lambda_ = lambda;
  seg_opts.enforce_consistent_boundaries_ = false;
  seg_opts.enforce_regionedge_ = false;
  seg_opts.base_filename_ = base_filename;

  if (app.is_set("-reduce-pairs"))
    seg_opts.reduce_edge_pairs_ = true;

  std::string constraint_string = app.getParam("-boundary-constraints");
  if (constraint_string == "tight") {
    seg_opts.enforce_regionedge_ = true;
  }
  else if (constraint_string == "simple") {
    seg_opts.enforce_consistent_boundaries_ = true;
  }
  else if (constraint_string == "extra") {
    seg_opts.enforce_regionedge_ = true;
    seg_opts.enforce_consistent_boundaries_ = true;
  }
  else if (constraint_string != "off") {
    USER_ERROR << "unknown boundary constraint mode \"" << constraint_string << "\"" << std::endl
	       << " choose from (tight | simple | extra | off)" << std::endl
	       << "  Exiting." << std::endl;
    exit(1);
  }

  if (gamma == 0.0) {
    
    lp_segment_lenreg(data_term, seg_opts, energy_offset, segmentation, &seeds);
  }
  else {
    seg_opts.gamma_ = gamma;
    seg_opts.bruckstein_ = app.is_set("-bruckstein");
    seg_opts.prevent_crossings_ = !app.is_set("-ignore-crossings");
    seg_opts.enforce_consistent_points_ = app.is_set("-no-touching-regions");
    seg_opts.light_constraints_ = app.is_set("-light-constraints");
    seg_opts.griddim_xDim_ = xDim;
    seg_opts.griddim_yDim_ = yDim;
    seg_opts.solver_ = app.getParam("-solver");

    int griddim = convert<int>(app.getParam("-griddim"));
    if (griddim > 0) {
      seg_opts.griddim_xDim_ = griddim;
      seg_opts.griddim_yDim_ = griddim;
    }
    seg_opts.debug_svg_ = app.is_set("-debug-svg");
    
    lp_segment_curvreg(data_term, seg_opts, energy_offset, segmentation ,&seeds);
  }

  tEndComputation = std::clock();

  std::cerr << "computation took " << diff_seconds(tEndComputation, tStartComputation) << " seconds." << std::endl;

  segmentation.savePGM(app.getParam("-o")+".out.ppm",255);
}
