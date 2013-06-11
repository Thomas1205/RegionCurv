/**** written by Petter Strandmark 2011 ****/

#ifndef ImageIntegrator_HH
#define ImageIntegrator_HH

#include <vector>
#include <string>

#include "mesh2D.hh"
#include "tensor.hh"
#include "matrix.hh"


class ImageIntegrator
{
  friend class SegmentationCurve;

public:
  ImageIntegrator(const Math2D::Matrix<float>& data_term);

  double integral(const std::vector<Mesh2DPoint>& coord) const;

private:
  double M(double x, int y) const;

  double fg_energy_line_pixel(double x1, double y1, double x2, double y2) const;
  double fg_energy_line(double x1, double y1, double x2, double y2) const;

  Math2D::Matrix<double> data_term_integrated_x;
  //Math2D::Matrix<double> data_term_integrated_y;
};



class SegmentationCurve
{
public:
  SegmentationCurve(const std::vector<Mesh2DPoint>& newcoord, ImageIntegrator& integrator_in, 
                    double sign, double lambda, double gamma, double curv_power, bool bruckstein);
  void reverse();

  bool self_intersect() const;
  void fix_self_intersect() ;
  void refine(bool verbose=false);

  double fg_energy();
  double smooth_energy();
  double energy();
  double energy2();

  double dEdx(int i, double h=1e-4);
  double dEdy(int i, double h=1e-4);

  double energy_single(int i, double x, double y);

  bool inside(size_t x, size_t y);

  static void start_svg(std::ofstream& of, const Math2D::Matrix<float>& image);
  static void end_svg(std::ofstream& of);
  void draw(std::ofstream& of);

private:

  bool step_cd();
  bool step_grad();

  ImageIntegrator* integrator;

  std::vector<Mesh2DPoint> coord;
  std::vector<Mesh2DPoint> original_coord;

  double lambda;
  double gamma; 
  double curv_power; 
  bool bruckstein;
};


#endif
