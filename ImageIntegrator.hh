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

	Math2D::Matrix<float> data_term_integrated_x;
  Math2D::Matrix<float> data_term_integrated_y;
};



class SegmentationCurve
{
public:
  SegmentationCurve(const std::vector<Mesh2DPoint>& newcoord, const ImageIntegrator& integrator_in, 
                    double lambda, double gamma, double curv_power, bool bruckstein);
	void reverse();

	void refine();

	double fg_energy();
	double smooth_energy();
	double energy();

	double dEdx(int i, double h=1e-4);
	double dEdy(int i, double h=1e-4);

	double energy_single(int i, double x, double y);

	void draw(std::string svgfile, const Math2D::Matrix<float>& image);

private:

	bool step_cd();
	bool step_grad();

  const ImageIntegrator& integrator;

	std::vector<Mesh2DPoint> coord;
	std::vector<Mesh2DPoint> original_coord;

	double lambda;
	double gamma; 
	double curv_power; 
	bool bruckstein;
};


#endif