/**** written by Petter Strandmark as an employee of Lund University, Sweden, 2010 ****/

#include <stdexcept>

#include "qpbo_segmentation.h"

#include "lp_segmentation.hh"
#include "segmentation_common.hh"
#include "mesh2D.hh"

#include "timing.hh"

#include "Petter-Color.hh"

#include "QPBO.h" 
#include "HOCR.h"

void err_function(char * err)
{
  //std::cerr << std::endl << Petter::RED << err << Petter::NORMAL << std::endl;
  throw std::runtime_error(err);
}


//solves a segmentation problem with length and curvature regularity via an LP
//@param lambda: the weight for the length term
//@param beta  : the weight for the curvature term  
double qpbo_segment_curvreg(const Math2D::Matrix<float>& data_term, const LPSegOptions& options, double energy_offset, Math2D::Matrix<uint>& segmentation)
{
  using namespace std;
  using namespace Petter;

  //Open the log file (append mode)
  string logfile_name = options.base_filename_ + ".qpbolog";
  ofstream logfile(logfile_name.c_str(), ios::app);
  logfile << options.lambda_ << " " << options.gamma_ << " ";


  const int svg_face_limit = 20000;

  const Math2D::Matrix<int>* fixed_labels = 0;

  //Get dimensions
  uint xDim = uint( data_term.xDim() );
  uint yDim = uint( data_term.yDim() );

  double lambda = options.lambda_;
  double gamma = options.gamma_;
  int neighborhood = options.neighborhood_;

  uint nAreasPerPixel = MAX_UINT;

  Mesh2D mesh;  
  if (options.gridtype_ == options.Square) {
    if (options.adaptive_mesh_n_ < 0) {
      double xfac = double(xDim) / double(options.griddim_xDim_); //Enlargement factor to get a griddim x griddim mesh
      double yfac = double(yDim) / double(options.griddim_yDim_); //Enlargement factor to get a griddim x griddim mesh
      generate_mesh( options.griddim_xDim_, options.griddim_yDim_, neighborhood,mesh);
      mesh.enlarge(xfac,yfac);
    }
    else {
      //Adaptive mesh
      generate_adaptive_mesh(data_term, mesh, neighborhood, options.adaptive_mesh_n_);
    }
  }
  else {
    if (options.adaptive_mesh_n_ < 0) {
      double xfac = double(xDim) / double(options.griddim_xDim_);
      double yfac = double(yDim) / double(options.griddim_yDim_);
      generate_hexagonal_mesh( xDim, yDim, 0.5*(xfac+yfac), neighborhood,mesh);
    }
    else {
      //Adaptive mesh
      generate_adaptive_hexagonal_mesh(data_term, mesh, neighborhood, options.adaptive_mesh_n_);
    }
  }

  std::vector<Mesh2DEdgePair> edge_pairs;
  statusTry("Generating edge pairs...");
  mesh.generate_edge_pair_list(edge_pairs);
  statusOK();

  std::cerr << edge_pairs.size() << " edge pairs." << std::endl;

  statusTry("Drawing mesh...");
  if (mesh.nFaces() > svg_face_limit) {
    statusFailed();
  }
  else {
    mesh.draw(options.base_filename_ + ".mesh.svg");
    statusOK();
  }

  statusTry("Calculating required QPBO size...");
  int required_nodes = mesh.nFaces();
  int required_edges = 0;
  for (uint j=0; j < edge_pairs.size(); j++) {

    uint first = edge_pairs[j].first_edge_idx_;
    uint second = edge_pairs[j].second_edge_idx_;

    if (mesh.adjacent_faces(first).size() != 2 || mesh.adjacent_faces(second).size() != 2)
    {
      //One of the edges are at the border of the image
      //No cost for this pair
      continue;
    }

    int weight = int(gamma  * curv_weight(mesh,edge_pairs[j],2,options.bruckstein_));

    if (weight == 0) {
      //This will not be requiring any extra nodes
      continue;
    }	

    uint first1 = mesh.adjacent_faces(first)[0];
    uint first2 = mesh.adjacent_faces(first)[1];
    uint second1 = mesh.adjacent_faces(second)[0];
    uint second2 = mesh.adjacent_faces(second)[1];

    //Vector of the first edge 
    double from1_x = mesh.point(mesh.edge(first).from_idx_).x_;
    double from1_y = mesh.point(mesh.edge(first).from_idx_).y_;
    double to1_x = mesh.point(mesh.edge(first).to_idx_).x_;
    double to1_y = mesh.point(mesh.edge(first).to_idx_).y_;

    //Vector of the second edge 
    double from2_x = mesh.point(mesh.edge(second).from_idx_).x_;
    double from2_y = mesh.point(mesh.edge(second).from_idx_).y_;
    double to2_x = mesh.point(mesh.edge(second).to_idx_).x_;
    double to2_y = mesh.point(mesh.edge(second).to_idx_).y_;

    //Make sure the second edge continues the first
    if (mesh.edge(first).from_idx_ == edge_pairs[j].common_point_idx_) {
      //Flip first edge
      double tmp = to1_x;
      to1_x = from1_x;
      from1_x = tmp;
      tmp = to1_y;
      to1_y = from1_y;
      from1_y = tmp;
    }
    if (mesh.edge(second).to_idx_ == edge_pairs[j].common_point_idx_) {
      //Flip second edge
      double tmp = to2_x;
      to2_x = from2_x;
      from2_x = tmp;
      tmp = to2_y;
      to2_y = from2_y;
      from2_y = tmp;
    }


    //cerr << endl << "Edge (" << from_x << "," << from_y << ") --> (" << to_x << "," << to_y << ")   faces " << first1 << " and " << first2 << endl;
    double ux = to1_x - from1_x;
    double uy = to1_y - from1_y;
    //Is first1 on the left or right side of the edge?
    Mesh2DPoint p = mesh.face_center(first1);
    //cerr << "First1 center = (" << p.x_ << "," << p.y_ << ") " << endl;
    double vx = p.x_  - from1_x;
    double vy = p.y_  - from1_y;
    double z = ux*vy - vx*uy;
    if (z > 0) {
      //This region is on the right side
      //Switch
      uint tmp = first1;
      first1 = first2;
      first2 = tmp;
      //cerr << "First1 on right side" << endl;
    }

    ux = to2_x - from2_x;
    uy = to2_y - from2_y;
    //Is second1 on the left or right side of the edge?
    p = mesh.face_center(second1);
    //cerr << "Second1 center = (" << p.x_ << "," << p.y_ << ") " << endl;
    vx = p.x_  - from2_x;
    vy = p.y_  - from2_y;
    z = ux*vy - vx*uy;
    if (z > 0) {
      //This region is on the right side
      //Switch
      uint tmp = second1;
      second1 = second2;
      second2 = tmp;
      //cerr << "Second1 on right side" << endl;
    }

    if (first1 == second1) {
      //Only three-clique needed
      required_edges += 3;
    }
    else if (first1 == second2) {
      //Not supposed to happen due to the switching above
      cerr << "first1 == second2" << endl;
      exit(1);
    }
    else if (first2 == second1) {
      //Not supposed to happen due to the switching above
      cerr << "first1 == second2" << endl;
      exit(1);
    }
    else if (first2 == second2) {
      required_edges += 3;
    }
    else {
      required_nodes += 10;
      required_edges += 24;
    }
  }

  for (uint j=0; j < mesh.nEdges(); j++) {
    if (mesh.adjacent_faces(j).size() != 2) {
      continue;
    }

    required_edges += 1;
  }

  statusOK();

  cerr << "Required QPBO nodes  : " << required_nodes << endl;
  cerr << "Required QPBO edges  : " << required_edges << endl;

  //
  // Create QPBO and HOCR optimizers
  //
  statusTry("Creating QPBO object...");
  QPBO<int> qpbo(required_nodes,required_edges, err_function);
  statusOK();
  statusTry("Creating HOCR object...");
  HOCR<int,4,QPBO<int> > hocr(qpbo);
  statusOK();
  statusTry("Adding nodes...");
  uint nFaces = mesh.nFaces();
  hocr.AddNode(nFaces);
  statusOK();

  //
  // Data terms for all regions
  //
  statusTry("Data weights... ");
  //for (uint y=0; y < yDim; y++) {
  //for (uint x=0; x < xDim; x++) {

  //	uint base = (y*xDim+x)*nAreasPerPixel;
  //	double cur_data = data_term(x,y); 

  //	for (uint i=0; i < nAreasPerPixel; i++) {
  //		double area = mesh.convex_area(base+i);
  //		double cost  = area * cur_data;
  //		double cost2 = calculate_data_term(base+i, mesh, data_term);
  //		if ( fabs(cost - cost2) > 1e-9) {
  //			cerr << "cost1=" << cost << " cost2=" << cost2 << endl;
  //		}
  //		hocr.AddUnaryTerm(base+i, 0, int(cost));
  //	}
  //}}

  float* raw_data = new float[mesh.nFaces()];
  for (int i=0; i<mesh.nFaces(); ++i) {
    double cost = calculate_data_term(i, mesh, data_term);
    raw_data[i] = float(cost);
    hocr.AddUnaryTerm(i, 0, int(cost));
  }
  statusOK();

  //
  // Length + curvature regularizer
  //

  QPBO<double>::EdgeId e;
  int nedges = 0;
  for (e=qpbo.GetNextEdgeId(-1); e>=0; e=qpbo.GetNextEdgeId(e))
  {
    nedges++;
  }
  cerr << "Number of QPBO nodes : " << qpbo.GetNodeNum() << endl;
  cerr << "Number of QPBO edges : " << nedges << endl;


  statusTry("Adding curvature weights...");


  int typical_weight = 0;
  double typical_curv_weight = 0;
  int nodes_added_4 = 0;
  int nodes_added_3 = 0;
  int arcs_added_4 = 0;
  int arcs_added_3 = 0;

  for (uint j=0; j < edge_pairs.size(); j++) {

    uint first = edge_pairs[j].first_edge_idx_;
    uint second = edge_pairs[j].second_edge_idx_;

    if (mesh.adjacent_faces(first).size() != 2 || mesh.adjacent_faces(second).size() != 2)
    {
      //One of the edges are at the border of the image
      //No cost for this pair
      continue;
    }

    int weight = int( gamma  * curv_weight(mesh,edge_pairs[j],2,options.bruckstein_));

    if (weight <= 0) {
      //This will not be requiring any treatment
      continue;
    }

    if (weight > 0 && typical_weight == 0) {
      typical_weight = weight;
      typical_curv_weight = curv_weight(mesh,edge_pairs[j]);
    }

    uint first1 = mesh.adjacent_faces(first)[0];
    uint first2 = mesh.adjacent_faces(first)[1];
    uint second1 = mesh.adjacent_faces(second)[0];
    uint second2 = mesh.adjacent_faces(second)[1];

    //Vector of the first edge 
    double from1_x = mesh.point(mesh.edge(first).from_idx_).x_;
    double from1_y = mesh.point(mesh.edge(first).from_idx_).y_;
    double to1_x = mesh.point(mesh.edge(first).to_idx_).x_;
    double to1_y = mesh.point(mesh.edge(first).to_idx_).y_;

    //Vector of the second edge 
    double from2_x = mesh.point(mesh.edge(second).from_idx_).x_;
    double from2_y = mesh.point(mesh.edge(second).from_idx_).y_;
    double to2_x = mesh.point(mesh.edge(second).to_idx_).x_;
    double to2_y = mesh.point(mesh.edge(second).to_idx_).y_;

    //Make sure the second edge continues the first
    if (mesh.edge(first).from_idx_ == edge_pairs[j].common_point_idx_) {
      //Flip first edge
      double tmp = to1_x;
      to1_x = from1_x;
      from1_x = tmp;
      tmp = to1_y;
      to1_y = from1_y;
      from1_y = tmp;
    }
    if (mesh.edge(second).to_idx_ == edge_pairs[j].common_point_idx_) {
      //Flip second edge
      double tmp = to2_x;
      to2_x = from2_x;
      from2_x = tmp;
      tmp = to2_y;
      to2_y = from2_y;
      from2_y = tmp;
    }


    //cerr << endl; 
    //cerr << "Edge (" << from1_x << "," << from1_y << ") --> (" << to1_x << "," << to1_y << ")" << endl;
    //cerr << "Edge (" << from2_x << "," << from2_y << ") --> (" << to2_x << "," << to2_y << ")" << endl;
    //cerr << "f1 = " << first1 << " f2 = " << first2 << " s1 = " << second1 << " s2 = " << second2 << endl; 

    //
    // Make sure that first1 and second1 are on the same side of the 
    // edge pair
    //

    //cerr << endl << "Edge (" << from_x << "," << from_y << ") --> (" << to_x << "," << to_y << ")   faces " << first1 << " and " << first2 << endl;
    double ux = to1_x - from1_x;
    double uy = to1_y - from1_y;
    //Is first1 on the left or right side of the edge?
    Mesh2DPoint p = mesh.face_center(first1);
    //cerr << "First1 center = (" << p.x_ << "," << p.y_ << ") " << endl;
    double vx = p.x_  - from1_x;
    double vy = p.y_  - from1_y;
    double z = ux*vy - vx*uy;
    if (z > 0) {
      //This region is on the right side
      //Switch
      uint tmp = first1;
      first1 = first2;
      first2 = tmp;
      //cerr << "First1 on right side" << endl;
    }

    ux = to2_x - from2_x;
    uy = to2_y - from2_y;
    //Is second1 on the left or right side of the edge?
    p = mesh.face_center(second1);
    //cerr << "Second1 center = (" << p.x_ << "," << p.y_ << ") " << endl;
    vx = p.x_  - from2_x;
    vy = p.y_  - from2_y;
    z = ux*vy - vx*uy;
    if (z > 0) {
      //This region is on the right side
      //Switch
      uint tmp = second1;
      second1 = second2;
      second2 = tmp;
      //cerr << "Second1 on right side" << endl;
    }

    //cerr << "f1 = " << first1 << " f2 = " << first2 << " s1 = " << second1 << " s2 = " << second2 << endl; 



    if (first1 == second1) {
      //Only three-clique needed
      int vals[8] = {0,0,0,weight, weight,0,0,0};
      int vars[3] = {first1, first2, second2};
      hocr.AddHigherTerm(3, vars, vals);
    }
    else if (first1 == second2) {
      //Not supposed to happen due to the switching above
      cerr << "first1 == second2" << endl;
      exit(1);
    }
    else if (first2 == second1) {
      //Not supposed to happen due to the switching above
      cerr << "first1 == second2" << endl;
      exit(1);
    }
    else if (first2 == second2) {
      //Only three-clique needed
      int vals[8] = {0,0,0,weight, weight,0,0,0};
      int vars[3] = {first2, first1, second1};

      int n = qpbo.GetNodeNum();
      int narcs = 0;
      if (arcs_added_3 == 0) {
        for (e=qpbo.GetNextEdgeId(-1); e>=0; e=qpbo.GetNextEdgeId(e)){
          narcs++;
        }	
      }

      hocr.AddHigherTerm(3, vars, vals);

      if (nodes_added_3 == 0) {
        nodes_added_3 = qpbo.GetNodeNum() - n;
      }

      if (arcs_added_3 == 0) {
        int narcs2 = 0;
        for (e=qpbo.GetNextEdgeId(-1); e>=0; e=qpbo.GetNextEdgeId(e)){
          narcs2++;
        }
        arcs_added_3 = narcs2-narcs;
      }
    }
    else {
      //We need a 4-clique to represent this energy
      //0011 -- cost
      //1100 -- cost
      //int vals[16] = {0,0,0,weight, 0,0,0,0, 0,0,0,0, weight,0,0,0}; //BUG?
      //0101 -- cost
      //1010 -- cost
      int vals[16] = {0,0,0,0, 0,weight,0,0, 0,0,weight,0, 0,0,0,0};
      int vars[4] = {first1, first2, second1, second2};

      int n = qpbo.GetNodeNum();
      int narcs = 0;
      if (arcs_added_4 == 0) {
        for (e=qpbo.GetNextEdgeId(-1); e>=0; e=qpbo.GetNextEdgeId(e)){
          narcs++;
        }	
      }

      hocr.AddHigherTerm(4, vars, vals);

      if (nodes_added_4 == 0) {
        nodes_added_4 = qpbo.GetNodeNum() - n;
      }
      if (arcs_added_4 == 0) {
        int narcs2 = 0;
        for (e=qpbo.GetNextEdgeId(-1); e>=0; e=qpbo.GetNextEdgeId(e)){
          narcs2++;
        }
        arcs_added_4 = narcs2-narcs;
      }
    }


    //cost[edge_pair_offset+2*j]   = weight;
    //cost[edge_pair_offset+2*j+1] = weight;
  }
  statusOK();

  statusTry("Adding length weights...");
  for (uint j=0; j < mesh.nEdges(); j++) {
    if (mesh.adjacent_faces(j).size() != 2) {
      continue;
    }

    int weight = int( lambda * mesh.edge_length(j) );
    uint r1 = mesh.adjacent_faces(j)[0];
    uint r2 = mesh.adjacent_faces(j)[1];
    qpbo.AddPairwiseTerm(r1,r2,0,weight,weight,0);
  }
  statusOK();

  nedges = 0;
  for (e=qpbo.GetNextEdgeId(-1); e>=0; e=qpbo.GetNextEdgeId(e))
  {
    nedges++;
  }




  cerr << "Number of QPBO nodes : " << qpbo.GetNodeNum() << endl;
  cerr << "Number of QPBO edges : " << nedges << endl;
  //cerr << "One positive pair cost: " << typical_weight << "  curvature weight : " << typical_curv_weight << endl;
//   cerr << "Nodes added for 4-clique : " << nodes_added_4 << " and " << arcs_added_4 << " arcs." << endl;
//   cerr << "Nodes added for 3-clique : " << nodes_added_3 << " and " << arcs_added_3 << " arcs." << endl;

  statusTry("Merging parallel edges....");
  qpbo.MergeParallelEdges();
  nedges = 0;
  for (e=qpbo.GetNextEdgeId(-1); e>=0; e=qpbo.GetNextEdgeId(e))
  {
    nedges++;
  }
  statusOK();
  cerr << "Number of QPBO edges : " << nedges << endl;

  //
  // Run optimizer
  //
  statusTry("Running QPBO...");

  qpbo.Solve();
  qpbo.ComputeWeakPersistencies();


  int unlabelled=0;
  for (uint i=0; i<mesh.nFaces(); ++i) {
    if (qpbo.GetLabel(i) < 0)
      unlabelled++;
  }

  int* labels = new int[mesh.nFaces()];
  for (uint i=0; i<mesh.nFaces(); ++i) {
    labels[i] = qpbo.GetLabel(i);
  }

  double time = statusOK();
  logfile << time << " ";

  cerr << "Unlabelled regions    : " << unlabelled << " (" << 100*double(unlabelled)/double(mesh.nFaces()) << "%)" << endl;
  logfile << 100*double(unlabelled)/double(mesh.nFaces()) << " ";
  bool no_regions_labeled = unlabelled == mesh.nFaces();

  unlabelled=0;
  for (int i=0; i<qpbo.GetNodeNum(); ++i) {
    if (qpbo.GetLabel(i) < 0)
      unlabelled++;
  }
  cerr << "Unlabelled QPBO nodes : " << unlabelled << " (" << 100*double(unlabelled)/double(qpbo.GetNodeNum()) << "%)" << endl;


  statusTry("Drawing SVG output...");
  
  if (mesh.nFaces() > svg_face_limit) {
    statusFailed();
  }
  else {
    mesh.draw_labels(options.base_filename_ + ".qpbo.out.svg",labels);
    statusOK();
  }

  // Only probe if we have unlabeled nodes.
  // Also, do not probe if *every* variable is unlableled -- it will not work
  // but potentially take a very long time.
  if (unlabelled > 0 /*&&  !no_regions_labeled*/ ) {
    int *mapping = new int[qpbo.GetNodeNum()];
    int *tmp_mapping = new int[qpbo.GetNodeNum()];
    for (int i = 0; i < qpbo.GetNodeNum(); i++) {
      mapping[i] = i * 2;
      tmp_mapping[i] = i * 2;
    }

    // Probe options
    QPBO<int>::ProbeOptions probeoptions;
    probeoptions.C = 200000;
    probeoptions.dilation = 1;
    probeoptions.weak_persistencies = 1;
    //options.iters = 1; //Only v.1.1


    statusTry("Running QPBO-P...");
    qpbo.Probe(mapping, probeoptions);
    qpbo.ComputeWeakPersistencies();
    /*for (int iter=1;iter<=5;++iter) {
    qpbo.Probe(tmp_mapping, options);
    qpbo.ComputeWeakPersistencies();
    qpbo.MergeMappings(nvars,mapping,tmp_mapping);
    }*/
    double time = statusOK();
    logfile << time << " ";

    // Read out entire labelling again (as weak persistencies may have changed)
    for (int nodeCount = 0; nodeCount < mesh.nFaces(); nodeCount++) {
      labels[nodeCount] = (int)qpbo.GetLabel(mapping[nodeCount]/2);
      if (labels[nodeCount] >= 0)
        labels[nodeCount] = (labels[nodeCount] + mapping[nodeCount]) % 2;
    }

    statusTry("Drawing SVG output...");
    if (mesh.nFaces() > svg_face_limit) {
      statusFailed();
    }
    else {
      mesh.draw_labels(options.base_filename_ + ".qpbop.svg",labels);
      statusOK();
    }

    //statusTry("Running QPBO-I...");
    //for (int iter=1;iter<=30;++iter) {
    //  qpbo.Improve();
    //}
    //statusOK();

    //// Read out entire labelling again (as weak persistencies may have changed)
    //for (int nodeCount = 0; nodeCount < nvars; nodeCount++) {
    //  labels[nodeCount] = (int)qpbo.GetLabel(mapping[nodeCount]/2);
    //  if (labels[nodeCount] >= 0)
    //    labels[nodeCount] = (labels[nodeCount] + mapping[nodeCount]) % 2;
    //}

    //delete[] mapping;
    //delete[] tmp_mapping;  

    //statusTry("Drawing SVG output...");
    //if (mesh.nFaces() > svg_face_limit) {
    //  statusFailed();
    //}
    //else {
    //  mesh.draw_labels(options.base_filename_ + ".qpboi.svg",labels);
    //  statusOK();
    //}
  }
  else {
    logfile << "-1 ";
  }


  statusTry("Building output...");

  //Calculate unlabelled pixels
  unlabelled = 0;
  for (uint i=0; i<mesh.nFaces(); ++i) {
    if (labels[i] < 0)
      unlabelled++;
  }
  //Reset segmentation
  Math2D::Matrix<double> output(xDim,yDim,0);
  for (uint i=0;i<mesh.nFaces();++i) {
    add_grid_output(i,labels[i],mesh,output);	
  }
  for (uint y=0; y < yDim; y++) {
    for (uint x=0; x < xDim; x++) {
      segmentation(x,y) = uint(output(x,y)*255.0);
  }}
  statusOK();

  delete[] raw_data;

  cerr << "Unlabelled regions    : " << unlabelled << " (" << 100*double(unlabelled)/double(mesh.nFaces()) << "%)" << endl;
  logfile << 100*double(unlabelled)/double(mesh.nFaces()) << " ";

  unlabelled=0;
  for (int i=0; i<qpbo.GetNodeNum(); ++i) {
    if (qpbo.GetLabel(i) < 0)
      unlabelled++;
  }
  cerr << "Unlabelled QPBO nodes : " << unlabelled << " (" << 100*double(unlabelled)/double(qpbo.GetNodeNum()) << "%)" << endl;


  delete[] labels;

  //Close logfile
  logfile << endl;

  return 0;
}


