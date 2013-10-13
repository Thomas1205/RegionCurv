/******* main program written by Thomas Schoenemann as an employee of the University of Pisa, Italy, 2011 *******/

#include "application.hh"
#include "colorimage.hh"
#include "curv_denoising.hh"

int main(int argc, char** argv) {

  if (argc == 1 || strings_equal(argv[1],"-h")) {

    std::cerr << "USAGE: " << argv[0] << std::endl
	      << "  -i <pgm or ppm> : filename of input image (to be denoised)" << std::endl
	      << "  -gamma <double> : curvature weight" << std::endl
	      << "  -lambda <double> : length weight" << std::endl
	      << "  -o <filename> : name of the (denoised) output image" << std::endl
	      << "  -n ( 4 | 8 | 16) : neighborhood size" << std::endl 
	      << "  -solver (clp | gurobi | xpress | cplex) : lp-solver used for curvature denoising" << std::endl
	      << "  [-b <uint>] : number of bins to approximate level lines (default: 1)" << std::endl
	      << "  [-hex-mesh] : use hexagonal meshes (in curvature mode)" << std::endl
	      << "  [-reduce-pairs] : disregard edge pairs with high curvature (in curvature mode, saves a little memory)" << std::endl 
	      << "  [-adaptive <uint> ] : number of regions in adaptive mode for curvature" << std::endl
	      << std::endl;
    exit(0);
  }

  const int nParams = 14;
  ParamDescr  params[nParams] = {{"-i",mandInFilename,0,""},{"-o",mandOutFilename,0,""},{"-n",optWithValue,1,"8"},
				 {"-gamma",mandWithValue,0,""}, {"-enforce-boundary-consistency",flag,0,""},
				 {"-enforce-region-edge-con",flag,0,""},{"-b",optWithValue,1,"1"},{"-bruckstein",flag,0,""},
				 {"-light-constraints",flag,0,""},{"-solver",optWithValue,1,"clp"},{"-lambda",optWithValue,1,"0.0"},
				 {"-hex-mesh",flag,0,""},{"-reduce-pairs",flag,0,""},{"-adaptive",optWithValue,0,""}};

  Application app(argc,argv,params,nParams);

  Math3D::NamedColorImage<float> input_image(app.getParam("-i"),MAKENAME(input_image));
  Math3D::NamedColorImage<double> out(MAKENAME(out));
  out.set_max_intensity(input_image.max_intensity());

  double gamma = convert<double>(app.getParam("-gamma"));

  LPSegOptions seg_opts;
  seg_opts.neighborhood_ = convert<uint>(app.getParam("-n"));
  seg_opts.lambda_ = convert<double>(app.getParam("-lambda"));
  seg_opts.gamma_ = gamma;
  seg_opts.bruckstein_ = app.is_set("-bruckstein");
  seg_opts.enforce_consistent_boundaries_ = app.is_set("-enforce-boundary-consistency");
  seg_opts.enforce_regionedge_ = app.is_set("-enforce-region-edge-con");
  seg_opts.griddim_xDim_ = input_image.xDim();
  seg_opts.griddim_yDim_ = input_image.yDim();
  seg_opts.light_constraints_ = app.is_set("-light-constraints");
  seg_opts.solver_ = app.getParam("-solver");
  if (app.is_set("-reduce-pairs"))
    seg_opts.reduce_edge_pairs_ = true;
  if (app.is_set("-hex-mesh"))
    seg_opts.gridtype_ = LPSegOptions::Hex;
  if (app.is_set("-adaptive"))
    seg_opts.adaptive_mesh_n_ = convert<int>(app.getParam("-adaptive"));
  
  curv_denoise(input_image, seg_opts, out, convert<uint>(app.getParam("-b")));    

  out.savePPM(app.getParam("-o"));
}
