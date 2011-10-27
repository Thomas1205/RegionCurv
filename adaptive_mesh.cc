/**** written by Petter Strandmark as an employee of Lund University, Sweden, summer 2010 ***/

#include <queue>
#include <string>
#include <fstream>

#include "gpcpetter.hh"

#include "mesh2D.hh"
#include "segmentation_common.hh"
#include "svg.hh"

#include "Petter-Color.hh"

//
// Basic element when generating adaptive square meshes
//
struct Rect
{
  size_t x1,x2,y1,y2;

  Rect(size_t x1,size_t x2,size_t y1, size_t y2) 
  {
    this->x1 = x1;
    this->x2 = x2;
    this->y1 = y1;
    this->y2 = y2;
  }
  size_t size()
  {
    return (x2-x1)*(y2-y1);
  }

  float compute_score(const Math2D::Matrix<float>& data_term)
  {
    float mean = 0;
    for (size_t x=x1; x<x2; ++x) {
      for (size_t y=y1; y<y2; ++y) {
        mean += data_term(x,y);
      }
    }
    mean /= size();

    float score = 0;
    for (size_t x=x1; x<x2; ++x) {
      for (size_t y=y1; y<y2; ++y) {
        float term = (data_term(x,y) - mean);
        score += term*term;
      }
    }

    return score;
  }
};

//
// Basic element when generating adaptive hexagonal meshes
//
struct Triangle
{
  Mesh2DPoint p1,p2,p3;

  Triangle(Mesh2DPoint p1, Mesh2DPoint p2, Mesh2DPoint p3) 
  {
    this->p1 = p1;
    this->p2 = p2;
    this->p3 = p3;
  }
  double size()
  {
    return triangle_area(p1,p2,p3);
  }

  int minx()
  {
    return int( floor(std::min(p1.x_,std::min(p2.x_,p3.x_))) + 0.5 );
  }
  int maxx()
  {
    return int( floor(std::max(p1.x_,std::max(p2.x_,p3.x_))) + 0.5 );
  }
  int miny()
  {
    return int( ceil(std::min(p1.y_,std::min(p2.y_,p3.y_))) + 0.5 );
  }
  int maxy()
  {
    return int( ceil(std::max(p1.y_,std::max(p2.y_,p3.y_))) + 0.5 );
  }
  float compute_score(const Math2D::Matrix<float>& data_term)
  {
    std::vector<Mesh2DPoint> points(3);
    points[0] = p1;
    points[1] = p2;
    points[2] = p3;

    //Polygon which represents a pixel
    std::vector<Mesh2DPoint> pixel(4);

    std::vector<double> areas;
    std::vector<double> values;
    for (int y=miny(); y < maxy(); ++y) {
      for (int x=minx(); x < maxx(); ++x) {
        if (x<0 || x>= ((int) data_term.xDim()) || y<0 || y>=((int) data_term.yDim())) continue;

        //Intersect the region with this pixel
        pixel[0].x_ = x;
        pixel[0].y_ = y;
        pixel[1].x_ = x+1;
        pixel[1].y_ = y;
        pixel[2].x_ = x+1;
        pixel[2].y_ = y+1;
        pixel[3].x_ = x;
        pixel[3].y_ = y+1;

        std::vector<Mesh2DPoint> intersection = poly_intersection(pixel, points);

        double area = 0;
        if (intersection.size() > 2) {
          for (uint i=1; i<intersection.size()-1; ++i) {
            area += triangle_area(intersection[0], intersection[i], intersection[i+1]);
          }
        }
        if (area > 1e-6) {
          areas.push_back(area);
          values.push_back(data_term(x,y));
        }
      }
    }

    double mean = 0;
    //Go through every pixel that possibly intersects
    //the region
    //double total = 0.0;
    double total_area = 0.0;
    for (uint i=0;i<areas.size(); ++i) {
      mean += areas[i]*values[i];
      total_area += areas[i];
    }
    mean /= total_area;

    double score = 0.0;
    for (uint i=0;i<areas.size(); ++i) {
      const double term = (values[i] - mean);
      score += areas[i]*term*term;
    }

    return float(score);
  }
};

//
// Comparator for Pairs.
// Pairs are compared using the second value -- the score.
//
template< typename FirstType, typename SecondType >
struct PairComparator {
  typedef const std::pair<FirstType, SecondType>& cref;
  bool operator()(cref p1, cref p2) const 
  {  
    return p1.second < p2.second;
  }
};


void add_rect_to_mesh(const Rect& rect, uint neighborhood, Mesh2D& mesh);
void add_triangle_to_mesh(const Triangle& tri, Mesh2D& mesh);

void generate_adaptive_mesh(const Math2D::Matrix<float>& data_term, Mesh2D& mesh, uint neighborhood, uint limit)
{
  using namespace std;
  using namespace Petter;

  statusTry("Generating adaptive mesh...");

  size_t xDim = data_term.xDim();
  size_t yDim = data_term.yDim();

  //Queue to determine which rectangle to split
  typedef pair<Rect,float> rectpair;
  priority_queue<rectpair, vector<rectpair>, PairComparator<Rect,float> > queue;


  //Add a rectangle to the priority queue
  Rect rect(0,xDim,0,yDim);
  queue.push( rectpair(rect,0.0f) );

  while (queue.size() < limit) {
    //Get the rectangle that needs to be divided the most
    const rectpair& rpair = queue.top();
    const Rect& rect = rpair.first;

    //Split into 4
    size_t xs = (rect.x2 + rect.x1)/2;
    size_t ys = (rect.y2 + rect.y1)/2;
    Rect rect1(rect.x1,xs,rect.y1,ys);
    Rect rect2(xs,rect.x2,rect.y1,ys);
    Rect rect3(rect.x1,xs,ys,rect.y2);
    Rect rect4(xs,rect.x2,ys,rect.y2);

    queue.pop();


    //Add to queue
    if (rect1.size() > 0) {
      queue.push( rectpair(rect1, rect1.compute_score(data_term)) );
    }
    if (rect2.size() > 0) {
      queue.push( rectpair(rect2, rect2.compute_score(data_term)) );
    }
    if (rect3.size() > 0) {
      queue.push( rectpair(rect3, rect3.compute_score(data_term)) );
    }
    if (rect4.size() > 0) {
      queue.push( rectpair(rect4, rect4.compute_score(data_term)) );
    }
  }

  //Draw as SVG
  /*priority_queue<rectpair, vector<rectpair>, PairComparator<Rect,float> > queue_copy = queue;
  string fg_style = "fill:#eeeeee;stroke-width:0.01;stroke:black";
  std::ofstream of("adaptive_mesh.svg");
  init_svg_file(of,xDim,yDim);
  while (!queue_copy.empty()) {
  rectpair& rpair = queue_copy.top();
  Rect& rect = rpair.first;
  vector<pair<double,double> > points;
  points.push_back( pair<double,double>(rect.x1,rect.y1) );
  points.push_back( pair<double,double>(rect.x2,rect.y1) );
  points.push_back( pair<double,double>(rect.x2,rect.y2) );
  points.push_back( pair<double,double>(rect.x1,rect.y2) );
  svg_draw_polygon(of,fg_style,points);
  queue_copy.pop();
  }
  finish_svg_file(of);*/

  statusOK();

  statusTry("Building connectivity...");
  //Generate the mesh
  //Mesh2D mesh;
  while (!queue.empty()) {
    const rectpair& rpair = queue.top();
    const Rect& rect = rpair.first;

    add_rect_to_mesh(rect,neighborhood,mesh);

    queue.pop();
  }
  statusOK();

  statusTry("Merging parallel edges...");
  mesh.merge_parallel_edges();
  statusOK();

  statusTry("Drawing...");
  if (mesh.nFaces() > 20000) {
    statusFailed();
  }
  else {
    mesh.draw("adaptive_mesh_sq_gen.svg", false);
    statusOK();
  }

  std::cerr << mesh.nPoints() << " points, " << mesh.nFaces() << " faces, " << mesh.nEdges() << " edges." << std::endl;
}




void add_rect_to_mesh(const Rect& rect, uint neighborhood, Mesh2D& mesh)
{
  using namespace std;
  using namespace Petter;

  uint nAreasPerPixel = 1;
  if (neighborhood == 4) {
    nAreasPerPixel = 1;
  }
  else if (neighborhood == 8) {
    nAreasPerPixel = 4;
  }
  else if (neighborhood == 16) {
    nAreasPerPixel = 32;
  }
  else {
    INTERNAL_ERROR << "invalid neighborhood \"" << neighborhood << "\". Exiting." << std::endl;
    exit(1);
  }

  double x1 = (double) rect.x1;
  double x2 = (double) rect.x2;
  double y1 = (double) rect.y1;
  double y2 = (double)rect.y2;
  double dx = double (x2-x1);
  double dy = double(y2-y1);

  Mesh2D smallmesh;
  generate_mesh(1,1, neighborhood, smallmesh, true);

  std::vector<uint> s2l(smallmesh.nPoints(), 0);
  for (uint i=0; i<smallmesh.nPoints(); ++i) {
    Mesh2DPoint p =  smallmesh.point(i);
    //Convert to large mesh coordinate system
    p.x_ = x1 + dx*p.x_;
    p.y_ = y1 + dy*p.y_;
    s2l[i] = mesh.find_or_add_point(p);
  }

  std::vector<uint> es2l(smallmesh.nEdges(), 0);
  for (uint i=0; i<smallmesh.nEdges(); ++i) {
    Mesh2DEdge e = smallmesh.edge(i);
    uint from = s2l[e.from_idx_];
    uint to   = s2l[e.to_idx_];
    es2l[i] = mesh.find_or_add_edge(to,from);
  }

  for (std::vector<Mesh2DFace>::const_iterator fi = smallmesh.face_start(); fi < smallmesh.face_end(); ++fi) {
    const Mesh2DFace& f = *fi;
    std::vector<uint> face;
    face.reserve(f.edge_idx_.size());
    for (uint j=0;j<f.edge_idx_.size();++j) {
      face.push_back( es2l[f.edge_idx_[j]] );
    }
    mesh.add_face(face);
  }
}


void add_triangle_to_mesh(const Triangle& tri, Mesh2D& mesh)
{
  using namespace std;

  std::vector<uint> p(7);
  std::vector<uint> l(12);

  p[0] = mesh.find_or_add_point(tri.p1);
  p[1] = mesh.find_or_add_point(tri.p2);
  p[2] = mesh.find_or_add_point(tri.p3);

  p[3] = mesh.find_or_add_point( Mesh2DPoint( (tri.p1.x_ + tri.p2.x_)/2, (tri.p1.y_ + tri.p2.y_)/2 ) );
  p[4] = mesh.find_or_add_point( Mesh2DPoint( (tri.p2.x_ + tri.p3.x_)/2, (tri.p2.y_ + tri.p3.y_)/2 ) );
  p[5] = mesh.find_or_add_point( Mesh2DPoint( (tri.p3.x_ + tri.p1.x_)/2, (tri.p3.y_ + tri.p1.y_)/2 ) );

  p[6] = mesh.find_or_add_point( Mesh2DPoint( (tri.p1.x_ + tri.p2.x_+ tri.p3.x_)/3, (tri.p1.y_ + tri.p2.y_ + tri.p3.y_)/3 ) );


  l[0] = mesh.find_or_add_edge(p[0],p[3]);
  l[1] = mesh.find_or_add_edge(p[3],p[1]);
  l[2] = mesh.find_or_add_edge(p[1],p[4]);
  l[3] = mesh.find_or_add_edge(p[4],p[2]);
  l[4] = mesh.find_or_add_edge(p[2],p[5]);
  l[5] = mesh.find_or_add_edge(p[5],p[0]);

  l[6]  = mesh.find_or_add_edge(p[6],p[0]);
  l[7]  = mesh.find_or_add_edge(p[6],p[3]);
  l[8]  = mesh.find_or_add_edge(p[6],p[1]);
  l[9]  = mesh.find_or_add_edge(p[6],p[4]);
  l[10] = mesh.find_or_add_edge(p[6],p[2]);
  l[11] = mesh.find_or_add_edge(p[6],p[5]);

  std::vector<uint> face;

  face.clear();
  face.push_back(l[6]);
  face.push_back(l[0]);
  face.push_back(l[7]);
  mesh.add_face(face);

  face.clear();
  face.push_back(l[7]);
  face.push_back(l[1]);
  face.push_back(l[8]);
  mesh.add_face(face);

  face.clear();
  face.push_back(l[8]);
  face.push_back(l[2]);
  face.push_back(l[9]);
  mesh.add_face(face);

  face.clear();
  face.push_back(l[9]);
  face.push_back(l[3]);
  face.push_back(l[10]);
  mesh.add_face(face);

  face.clear();
  face.push_back(l[10]);
  face.push_back(l[4]);
  face.push_back(l[11]);
  mesh.add_face(face);

  face.clear();
  face.push_back(l[11]);
  face.push_back(l[5]);
  face.push_back(l[6]);
  mesh.add_face(face);
}


void generate_adaptive_hexagonal_mesh(const Math2D::Matrix<float>& data_term, Mesh2D& mesh, uint neighborhood, uint limit)
{
  using namespace std;
  using namespace Petter;
  statusTry("Generating triangular mesh...");

  Mesh2DPoint p1(-double(data_term.xDim())/2, double(data_term.yDim()));
  Mesh2DPoint p2(3*double(data_term.xDim())/2, double(data_term.yDim()));
  Mesh2DPoint p3(double(data_term.xDim())/2, -double(data_term.yDim())/2);

  //Queue to determine which rectangle to split
  typedef pair<Triangle,float> tripair;
  std::priority_queue<tripair, vector<tripair>, PairComparator<Triangle,float> > queue;


  //Add a first triangle to the priority queue
  Triangle tri(p1,p2,p3);
  queue.push( tripair(tri,0.0f) );

  while (queue.size() < limit) {
    //Get the rectangle that needs to be divided the most
    const tripair& rpair = queue.top();
    const Triangle& tri = rpair.first;

    //Split into 4
    Triangle tri1( tri.p1, 
      Mesh2DPoint( (tri.p1.x_ + tri.p2.x_)/2 , (tri.p1.y_ + tri.p2.y_)/2 ),
      Mesh2DPoint( (tri.p1.x_ + tri.p3.x_)/2 , (tri.p1.y_ + tri.p3.y_)/2 ) );
    Triangle tri2( tri.p2, 
      Mesh2DPoint( (tri.p2.x_ + tri.p3.x_)/2 , (tri.p2.y_ + tri.p3.y_)/2 ),
      Mesh2DPoint( (tri.p2.x_ + tri.p1.x_)/2 , (tri.p2.y_ + tri.p1.y_)/2 ) );
    Triangle tri3( tri.p3, 
      Mesh2DPoint( (tri.p3.x_ + tri.p1.x_)/2 , (tri.p3.y_ + tri.p1.y_)/2 ),
      Mesh2DPoint( (tri.p3.x_ + tri.p2.x_)/2 , (tri.p3.y_ + tri.p2.y_)/2 ) );
    Triangle tri4( Mesh2DPoint( (tri.p1.x_ + tri.p2.x_)/2 , (tri.p1.y_ + tri.p2.y_)/2 ), 
      Mesh2DPoint( (tri.p2.x_ + tri.p3.x_)/2 , (tri.p2.y_ + tri.p3.y_)/2 ),
      Mesh2DPoint( (tri.p3.x_ + tri.p1.x_)/2 , (tri.p3.y_ + tri.p1.y_)/2 ) );

    queue.pop();


    //Add to queue
    queue.push( tripair(tri1, tri1.compute_score(data_term)) );
    queue.push( tripair(tri2, tri2.compute_score(data_term)) );
    queue.push( tripair(tri3, tri3.compute_score(data_term)) );
    queue.push( tripair(tri4, tri4.compute_score(data_term)) );
  }

  statusOK();

  statusTry("Building connectivity...");
  //Generate the mesh
  //Mesh2D mesh;
  while (!queue.empty()) {
    const tripair& rpair = queue.top();
    const Triangle& tri = rpair.first;

    add_triangle_to_mesh(tri,mesh);

    queue.pop();
  }
  statusOK();

  statusTry("Merging parallel edges...");
  mesh.merge_parallel_edges();
  statusOK();

  statusTry("Drawing...");
  if (mesh.nFaces() > 20000) {
    statusFailed();
  }
  else {
    mesh.draw("adaptive_mesh_hex_gen.svg", false);
    statusOK();
  }

  std::cerr << mesh.nPoints() << " points, " << mesh.nFaces() << " faces, " << mesh.nEdges() << " edges." << std::endl;
}
