/*** first version written by Thomas Schoenemann as a private person without employment, September 2009 ***/
/*** continued by Thomas Schoenemann as an employee of Lund University, Sweden, 2010 - 2011 ***/
/***  with additional support by Petter Strandmark ***/

#include "application.hh"
#include "grayimage.hh"
#include "colorimage.hh"
#include "conversion.hh"
#include "draw_segmentation.hh"
#include "lp_segmentation.hh"
#include "extended_lp_segmentation.hh"
#include "timing.hh"
#include "sampling.hh"
#include "lp_segmenter.hh"

#ifdef HAS_QPBO
//by Petter Strandmark
#include "qpbo_segmentation.hh"
#endif

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
      << "  -i <pgm or ppm> : filename of input image (to be segmented)" << std::endl
      << "  -lambda <double> : length weight" << std::endl
      << "  -gamma <double> : curvature weight" << std::endl
      << "  -o <filename> : name of the output segmentation" << std::endl
      << "  -method (lp | factor-lp) " << std::endl
      << "  -regions <uint>" << std::endl
      << "  -mu0 <double> -mu1 <double> : mean values for background and foreground, respectively" << std::endl
      << "  -griddim <uint> : dimension of the mesh, default is the image size" << std::endl
      << "     if you want a non-square mesh different from the image size, use -griddimx <uint> and -griddimy <uint>" << std::endl
      << "  -boundary-constraints (tight | simple | extra | off) : constraints to ensure consistency of regions and boundaries." 
      << "     default is tight (= Strandmark&Kahl 11), extra unites simple and tight " << std::endl 
      << " [-n (4|8|16)] : size of neighborhood, default 8" << std::endl
      << " [-ignore-crossings] : allow crossings of line pairs, e.g. when three regions meet in a point" << std::endl
      << " [-light-constraints] : take only half of the optional constraints" << std::endl
      << " [-hex-mesh] : use hexagonal mesh instead of squared mesh" << std::endl
      << " [-adaptive (uint)] : use adaptive meshing" << std::endl
      << " [-debug-svg]: draw SVG files for debugging" << std::endl
      << " [-reduce-pairs] : same some memory and run-time by not considering pairs with very high curvature" << std::endl
      << " -solver ( clp | gurobi | cplex | xpress | own-conv ) : default clp" << std::endl
      << " -mode ( standard | msed | qpbo | icm )" << std::endl
      << std::endl;

    exit(0);
  }

  ParamDescr  params[] = {{"-i",mandInFilename,0,""},{"-lambda",optWithValue,1,"1.0"},
  {"-gamma",optWithValue,1,"1.0"},{"-o",mandOutFilename,0,""},{"-n",optWithValue,1,"8"},
  {"-boundary-constraints",optWithValue,1,"tight"},
  {"-method",optWithValue,1,"lp"},{"-bruckstein",flag,0,""},
  {"-adaptive",optWithValue,0,""},{"-hex-mesh",flag,0,""},{"-regions",optWithValue,1,"2"},
  {"-light-constraints",flag,0,""},{"-debug-svg",flag,0,""},{"-mu0",optWithValue,1,"-1"},
  {"-mu1",optWithValue,1,"-1"},{"-griddim",optWithValue,1,"-1"},{"-griddimx",optWithValue,1,"-1"},
  {"-griddimy",optWithValue,1,"-1"},{"-solver",optWithValue,1,"clp"},
  {"-ignore-crossings",flag,0,""},{"-no-touching-regions",flag,0,""},{"-convex",flag,0,""},
  {"-reduce-pairs",flag,0,""},{"-mode",optWithValue,1,"standard"},
  {"-refine",flag,0,""}};

  const int nParams = sizeof(params)/sizeof(ParamDescr);

  Application app(argc,argv,params,nParams);

  std::string base_filename = app.getParam("-o");
  //check_filename(base_filename + ".final.svg");
  //check_filename(base_filename + ".lp.svg");
  //check_filename(base_filename + ".lp_simple.svg");
  check_filename(base_filename + ".out.ppm");
  check_filename(base_filename + ".seg.pgm");
  //check_filename(base_filename + ".frac.pgm");

  Math3D::NamedColorImage<float> color_image(app.getParam("-i"),MAKENAME(color_image));

  uint xDim = uint( color_image.xDim() );
  uint yDim = uint( color_image.yDim() );
  uint zDim = uint( color_image.zDim() );

  Math2D::NamedGrayImage<float> gray_image(xDim,yDim,color_image.max_intensity(),MAKENAME(gray_image));
  for (uint y=0; y < yDim; y++) {
    for (uint x=0; x < xDim; x++) {

      double sum = 0.0;
      for (uint z=0; z < zDim; z++)
        sum += color_image(x,y,z);

      gray_image(x,y) = float(sum / zDim);
    }
  }

  uint nRegions = convert<uint>(app.getParam("-regions"));
  uint neighborhood = convert<uint>(app.getParam("-n"));

  std::string method_string = app.getParam("-method");
  std::string mode_string = app.getParam("-mode");

  Math2D::NamedMatrix<float> data_term(xDim,yDim,MAKENAME(data_term));
  Math2D::NamedMatrix<uint> segmentation(xDim,yDim,0,MAKENAME(segmentation));

  float energy_offset = 0.0;
  float mu0 = convert<float>(app.getParam("-mu0"));
  float mu1 = convert<float>(app.getParam("-mu1"));
  if (mu0 < 0) {
    mu0 = gray_image.min();
  }
  if (mu1 < 0) {
    mu1 = gray_image.max();
  }

  for (uint y=0; y < yDim; y++) {
    for (uint x=0; x < xDim; x++) {
      float cur = gray_image(x,y);
      float data0 = (cur-mu0)*(cur-mu0);
      float data1 = (cur-mu1)*(cur-mu1);
      data_term(x,y) = data1 - data0;
      energy_offset += data0;
    }
  }

  Math3D::NamedTensor<float> multi_data_term(MAKENAME(multi_data_term));

  if (nRegions > 2) {

    multi_data_term.resize(xDim,yDim,nRegions);

    Math1D::Vector<double> mean(nRegions,0.0);

    Math1D::Vector<float> intensity(gray_image.size());
    for (uint i=0; i < gray_image.size(); i++)
      intensity[i] = gray_image.direct_access(i);
    std::sort(intensity.direct_access(),intensity.direct_access()+gray_image.size());

    //spread percentiles equally
    double share = 1.0 / nRegions;
    for (uint i=0; i < nRegions; i++) {
      mean[i] = intensity[std::size_t(gray_image.size() * (i+0.5)*share)];
    }

    for (uint y=0; y < yDim; y++) {
      for (uint x=0; x < xDim; x++) {
        float cur = gray_image(x,y);
        for (uint r=0; r < nRegions; r++) {
          float data = float((cur-mean[r])*(cur-mean[r]));
          multi_data_term(x,y,r) = data;
        }
      }
    }
  }

  double lambda = convert<double>(app.getParam("-lambda"));
  double gamma = convert<double>(app.getParam("-gamma"));

  std::clock_t tStartComputation, tEndComputation;
  tStartComputation = std::clock();

  LPSegOptions seg_opts;
  seg_opts.neighborhood_ = neighborhood;
  seg_opts.lambda_ = lambda;
  seg_opts.gamma_ = gamma;
  seg_opts.bruckstein_ = app.is_set("-bruckstein");
  seg_opts.enforce_consistent_boundaries_ = false;
  seg_opts.enforce_regionedge_ = false;
  seg_opts.enforce_consistent_points_ = app.is_set("-no-touching-regions");
  seg_opts.prevent_crossings_ = !app.is_set("-ignore-crossings");
  seg_opts.light_constraints_ = app.is_set("-light-constraints");
  seg_opts.griddim_xDim_ = xDim;
  seg_opts.griddim_yDim_ = yDim;
  seg_opts.solver_ = app.getParam("-solver");
  seg_opts.base_filename_ = base_filename;
  seg_opts.convex_prior_ =  app.is_set("-convex");
  seg_opts.refine_ =  app.is_set("-refine");

  seg_opts.debug_svg_ = app.is_set("-debug-svg");

  if (app.is_set("-reduce-pairs"))
    seg_opts.reduce_edge_pairs_ = true;
  if (app.is_set("-hex-mesh"))
    seg_opts.gridtype_ = LPSegOptions::Hex;
  if (app.is_set("-adaptive"))
    seg_opts.adaptive_mesh_n_ = convert<int>(app.getParam("-adaptive"));

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

  //The dimensions of the mesh may be different form the image
  int griddim =  convert<int>(app.getParam("-griddim"));
  int griddim_x =  convert<int>(app.getParam("-griddimx"));
  int griddim_y =  convert<int>(app.getParam("-griddimy"));
  if (griddim_x < 0)
    griddim_x = griddim;
  if (griddim_y < 0)
    griddim_y = griddim;  

  if (griddim_x > 0 && griddim_y > 0) {
    seg_opts.griddim_xDim_ = griddim_x;
    seg_opts.griddim_yDim_ = griddim_y;
    seg_opts.output_factor_ = 1; // Usually no need to enlarge output now
  }

  if (gamma == 0.0 && !seg_opts.convex_prior_) {

    if (nRegions == 2)
      lp_segment_lenreg(data_term, seg_opts, energy_offset, segmentation);
    else {

      LpSegmenter lp_segmenter(gray_image, seg_opts, nRegions, false);

      lp_segmenter.segment(1);

      segmentation = lp_segmenter.segmentation();
    }
  }
  else {

    if (method_string == "lp") {
      if (nRegions == 2) {
        if (mode_string == "icm")
          curv_icm(data_term, seg_opts, energy_offset, segmentation);
        else
          lp_segment_curvreg(data_term, seg_opts, energy_offset, segmentation);
      }
      else {

        if (mode_string == "icm") {
          multi_curv_icm(multi_data_term, seg_opts, segmentation);
        }
        else {
          LpSegmenter lp_segmenter(gray_image, seg_opts, nRegions, false);

          lp_segmenter.segment(1);

          segmentation = lp_segmenter.segmentation();
        }
      }
    }
    else {

      if (mode_string == "msd") {
        //clique_lp_segment_curvreg_minsum_diffusion(data_term, seg_opts, energy_offset, segmentation);
        factor_lp_segment_curvreg_minsum_diffusion_memsave(multi_data_term, seg_opts, energy_offset, segmentation);
      }
      else if (mode_string == "qpbo") {
#ifdef HAS_QPBO
        qpbo_segment_curvreg(data_term, seg_opts, energy_offset, segmentation);
#else
        USER_ERROR << "method QPBO is not available in this pacakage" << std::endl;
        exit(1);
#endif
      }
      else if (mode_string == "icm") {
        if (nRegions == 2)
          factor_curv_icm(data_term, seg_opts, energy_offset, segmentation);
        else
          factor_curv_icm(multi_data_term, seg_opts, segmentation);
      }
      else {
        if (nRegions == 2)
          factor_lp_segment_curvreg(data_term, seg_opts, energy_offset, segmentation);
        else {
          //clique_lp_segment_pottscurvreg(multi_data_term, seg_opts, segmentation);
          factor_lp_segment_pottscurvreg_layered(multi_data_term, seg_opts, segmentation);
        }
      }
    }
  }


  tEndComputation = std::clock();

  std::cerr << "computation took " << diff_seconds(tEndComputation, tStartComputation) << " seconds." << std::endl;

  segmentation.savePGM(base_filename + ".seg.pgm",255);

  Math3D::NamedColorImage<float> out_image(MAKENAME(out_image));
  make_color_image(color_image,out_image);  

  uint scale_fac = 255 / (nRegions - 1);

  Math2D::Matrix<float> true_seg(segmentation.xDim(),segmentation.yDim());
  for (uint i=0; i < true_seg.size(); i++) {

    true_seg.direct_access(i) = double(segmentation.direct_access(i)) / double(scale_fac) + 0.5;
  }

  Math2D::Matrix<float> fscaled_seg(xDim,yDim);
  downsample_matrix(true_seg, fscaled_seg);

  Math2D::Matrix<uint> scaled_seg(xDim,yDim);
  for (uint i=0; i < scaled_seg.size(); i++) {
    scaled_seg.direct_access(i) = uint (fscaled_seg.direct_access(i) + 0.5);
  }

  draw_segmentation(scaled_seg, out_image, 250.0, 250.0, 150.0);

  out_image.savePPM(base_filename + ".out.ppm");
}
