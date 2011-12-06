/**** written by Petter Strandmark 2011 ****/



#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <cmath>
using namespace std;

#include "ImageIntegrator.hh"
#include "svg.hh"
#include "curvature.hh"

#define ASSERT(cond) if (!(cond)) { cout << endl << "Error (line " << __LINE__ << " in " << __FILE__ << "): " << #cond << endl; \
                                    exit(128); }

ImageIntegrator::ImageIntegrator(const Math2D::Matrix<float>& data_term) :
	data_term_integrated_x(data_term.xDim(), data_term.yDim()),
	data_term_integrated_y(data_term.xDim(), data_term.yDim())
{
  
	size_t xDim = data_term_integrated_x.xDim();
	size_t yDim = data_term_integrated_y.yDim();

	for (size_t y=0;y<yDim;++y) {
	for (size_t x=0;x<xDim;++x) {
		data_term_integrated_x(x,y) = data_term(x,y);
		data_term_integrated_y(x,y) = data_term(x,y);
	}}
	
	// Integrate data term
	for (size_t y=0;y<yDim;++y) {
		for (size_t x=1;x<xDim;++x) {
			data_term_integrated_x(x,y) += data_term_integrated_x(x-1,y);
		}
	}

	for (size_t x=0;x<xDim;++x) {
		for (size_t y=1;y<yDim;++y) {
			data_term_integrated_y(x,y) += data_term_integrated_y(x,y-1);
		}
	}
}



//
// M is linear in x and constant in y over each pixel
//
double ImageIntegrator::M(double px, int y) const
{

	// Outside image region, the image is constant
	if (px >= data_term_integrated_x.xDim()) {
		px = double(data_term_integrated_x.xDim());
	}
	if (y < 0) {
		y = 0;
	}
	if (y >= data_term_integrated_x.yDim()) {
		y = data_term_integrated_x.yDim()-1;
	}
	
	int x = int(px);
	//int y = int(py);
	double frac = px - int(px);
	
	int x0 = x - 1;
	int x1 = x;
	double f0,f1;
	if (x0 >= 0) {
		f0 = data_term_integrated_x(x0,y);
	} else {
		f0 = 0;
	}
	f1 = data_term_integrated_x(x1,y);

	return (1-frac)*f0 + frac*f1;
}


double ImageIntegrator::fg_energy_line_pixel(double x1, double y1, double x2, double y2) const
{
	double dy = y2 - y1;
	int y = int( (y1+y2)/2 );
	double Mavg = (M(x1,y) + M(x2,y))/2;
	return dy * Mavg;
}


double ImageIntegrator::fg_energy_line(double x1, double y1, double x2, double y2) const
{

	// Split the line into segments for each pixel
	// x = (1-t)*x1 + t*x2
	// y = (1-t)*y1 + t*y2
	vector<double> ts;
	ts.push_back(0);
	if (abs(x2-x1) > 0) {
		for (int x=int(min(x1,x2)); x<=int(max(x1,x2)); ++x) {
			ts.push_back( (x-x1) / (x2-x1) );
		}
	}
	if (abs(y2-y1) > 0) {
		for (int y=int(min(y1,y2)); y<=int(max(y1,y2)); ++y) {
			ts.push_back( (y-y1) / (y2-y1) );
		}
	}
	ts.push_back(1);

	sort(ts.begin(), ts.end());

	// Compute the energy for each line segment separately
	double en = 0;
	
	size_t xDim = data_term_integrated_x.xDim();
	size_t yDim = data_term_integrated_y.yDim();

  
	for (size_t i=0; i < ts.size()-1; ++i) {
		if (ts[i] >= 0 && ts[i+1] <= 1) {
			double xa = (1-ts[i])*x1 + ts[i]*x2;
			double ya = (1-ts[i])*y1 + ts[i]*y2;
			double xb = (1-ts[i+1])*x1 + ts[i+1]*x2;
			double yb = (1-ts[i+1])*y1 + ts[i+1]*y2;
			//Project points to image domain
			if (xa<0) xa = 0;
			if (xa>xDim) xa=xDim;
			if (ya<0) ya = 0;
			if (ya>yDim) ya=yDim;
			if (xb<0) xb = 0;
			if (xb>xDim) xb=xDim;
			if (yb<0) yb = 0;
			if (yb>yDim) yb=yDim;
			if ( (xa-xb)*(xa-xb) + (ya-yb)*(ya-yb) > 0 ) {
				en += fg_energy_line_pixel(xa,ya,xb,yb);
			}
		}
	}

	return en;
}

double ImageIntegrator::integral(const std::vector<Mesh2DPoint>& coord) const
{
	double en = 0;
	for (uint i=0; i < coord.size(); i++) {
		const Mesh2DPoint& p1 = coord[i];
		const Mesh2DPoint& p2 = coord[ (i+1) % coord.size() ];
		en += fg_energy_line(p1.x_, p1.y_, p2.x_, p2.y_);
	}
	return en;
}








SegmentationCurve::SegmentationCurve(const std::vector<Mesh2DPoint>& newcoord, 
                                     const ImageIntegrator& integrator_in, 
                                     double lambda, 
                                     double gamma, 
                                     double curv_power, 
                                     bool bruckstein) :
	coord(newcoord),
	original_coord(newcoord),
  integrator(integrator_in)
{
	this->lambda = lambda;
	this->gamma = gamma;
	this->curv_power = curv_power;
	this->bruckstein = bruckstein;
}

void SegmentationCurve::reverse()
{
	std::reverse(coord.begin(), coord.end());
	std::reverse(original_coord.begin(), original_coord.end());
}

bool SegmentationCurve::step_cd()
{
	bool any_changed = false;

	for (int i=0; i < coord.size(); i++) {
		//cerr << "i=" << setw(4) << i << " " << endl;

		double x = coord.at(i).x_;
		double y = coord.at(i).y_;
		double E1 = energy_single(i,x,y);
		double dx = dEdx(i);
		double dy = dEdy(i);
		//double dx2 = dEdx_ad(i);
		//double dy2 = dEdy_ad(i);

		//cerr << "dx=" << setw(20) << dx << " dy=" << setw(20) << dy << " " << endl;
		//cerr << "dx=" << setw(20) << dx2 << " dy=" << setw(20) << dy2 << " " << endl;

		double nrm = sqrt(dx*dx + dy*dy);

		if (nrm>1) {
			dx /= nrm;
			dy /= nrm;
		}

		double E2 = E1;
		double alpha = 1.1;
		double tau = 0.05;

		while (true) {
			double xnew = x - alpha*tau*dx;
			double ynew = y - alpha*tau*dy;
			double Enew = energy_single(i,xnew,ynew);

			if (Enew >= E2) {
				break;
			}

			tau = alpha*tau;
			E2 = Enew;

			if (tau > 0.4) {
				break;
			}
		} 

		//cerr << setw(20) << E1 << " --> " << setw(20) << E2;

		if (E2 < E1) {

			double Ecomplete1 = energy();
			double Esingle1   = energy_single(i,x,y);

			coord.at(i).x_ = x - tau*dx;
			coord.at(i).y_ = y - tau*dy;	

			double Ecomplete2 = energy();
			double Esingle2   = energy_single(i,x - tau*dx,y - tau*dy);

			ASSERT( (Ecomplete1 - Ecomplete2)/abs(Ecomplete1) > -1e-5 );
			ASSERT( (Esingle1 - Esingle2)/abs(Esingle1) > 0 );

			//cerr << " changed";
			any_changed = true;
		}
		//cerr << endl;
	}

	return any_changed;
}

bool SegmentationCurve::step_grad()
{
	// Compute gradient
	vector<double> dx(coord.size(), 0);
	vector<double> dy(coord.size(), 0);
	double max_nrm = -1;
	for (int i=0; i < coord.size(); i++) {
		dx[i] = dEdx(i);
		dy[i] = dEdy(i);
		double nrm = sqrt(dx[i]*dx[i] + dy[i]*dy[i]);
		if (nrm > max_nrm) {
			max_nrm = nrm;
		}
	}

	for (int i=0; i < coord.size(); i++) {
		dx[i] /= max_nrm;
		dy[i] /= max_nrm;
	}
	


	double E1 = energy();
	double E2 = E1;
	double alpha = 1.1;
	double tau = 0.05;

	vector<Mesh2DPoint> coordold(coord);
	vector<Mesh2DPoint> coordgood(coord);

	while (true) {
		for (int i=0; i < coord.size(); i++) {
			coord[i].x_ = coordold[i].x_ - alpha*tau*dx[i];
			coord[i].y_ = coordold[i].y_ - alpha*tau*dy[i];
		}
		double Enew = energy();
		if (Enew >= E2) {
			break;
		}
		if (tau > 1.0) {
			break;
		}
		tau = tau*alpha;
		E2 = Enew;
		coordgood = coord;
	} 

	coord = coordgood;

	return E2 < E1;
}


void SegmentationCurve::refine()
{
	lambda = 0;

	for (int iter=1; iter <= 200; ++iter) {

		double Estart = energy();
		std::vector<Mesh2DPoint> oldcoord = coord;

		bool any_changed = false;
		any_changed = step_cd() || any_changed;
		//any_changed = step_grad() || any_changed;

		double Eend = energy();

		//cerr << "--------------------" << endl;
		cerr << "iter=" << setw(3) << iter << " E = " << Eend << " ";
		if (Estart == Eend) {
			cerr << "(no change) ";
		}
		else if (Estart < Eend)  {
			cerr << "(increase) ";
		}
		else if (Estart > Eend) {
			cerr << "(decrease) ";
		}

		// Check that the two ways of computing the energy give
		// the same result
		double E2=0;
		for (int i=0; i < coord.size(); i++) {
			double x = coord.at(i).x_;
			double y = coord.at(i).y_;
			E2 += energy_single(i,x,y);
		}
		cerr << "deltaE = " << abs((E2-Eend)/E2) << " ";
		cerr << endl;
		//cerr << "--------------------" << endl;

		if (!any_changed || Estart < Eend) {
			break;
			coord = oldcoord;
		}
	}
}



double SegmentationCurve::fg_energy() 
{
  return integrator.integral(coord);
}

double SegmentationCurve::smooth_energy() 
{
	double en = 0;
	for (int i=0; i < coord.size(); i++) {
		int iprev = i-1;
		if (iprev<0) iprev = coord.size()-1;
		Mesh2DPoint& p0 = coord[ iprev ];
		Mesh2DPoint& p1 = coord[  i ];
		Mesh2DPoint& p2 = coord[ (i+1) % coord.size() ];

		double dx = p1.x_-p0.x_;
		double dy = p1.y_-p0.y_;
		en += 0.5 * lambda * sqrt( dx*dx + dy*dy );
		dx = p1.x_-p2.x_;
		dy = p1.y_-p2.y_;
		en += 0.5 * lambda * sqrt( dx*dx + dy*dy );

		en += gamma * curv_weight(p0.x_,p0.y_, p1.x_,p1.y_, p2.x_,p2.y_, 2.0, bruckstein);
		//cerr << p0 << " " << p1 << " " << p2 << "  en=" << en << endl;
	}

	return en;
}

double SegmentationCurve::energy() 
{
	return fg_energy() + smooth_energy();
}


double SegmentationCurve::dEdx(int i, double h)
{
	double x = coord.at(i).x_;
	double y = coord.at(i).y_;
	return ( energy_single(i, x + h, y) - energy_single(i, x - h, y) ) / (2*h);
}

double SegmentationCurve::dEdy(int i, double h)
{
	double x = coord.at(i).x_;
	double y = coord.at(i).y_;
	return ( energy_single(i, x, y + h) - energy_single(i, x, y - h) ) / (2*h);
}


double SegmentationCurve::energy_single(int i, double x, double y)
{
	double en=0;

	int iprev2 = i-2;
	if (iprev2<0) iprev2 += coord.size();
	int iprev = i-1;
	if (iprev<0) iprev = coord.size()-1;
	
	//   p0---p1---(x,y)---p2---p3

	Mesh2DPoint& p0 = coord[ iprev2 ];
	Mesh2DPoint& p1 = coord[ iprev ];
	Mesh2DPoint& p2 = coord[ (i+1) % coord.size() ];
	Mesh2DPoint& p3 = coord[ (i+2) % coord.size() ];

	double dx = x-p1.x_;
	double dy = y-p1.y_;
	en += 0.5 * lambda * sqrt( dx*dx + dy*dy );
	dx = x-p2.x_;
	dy = y-p2.y_;
	en += 0.5 * lambda * sqrt( dx*dx + dy*dy );

	en += 0.33333333333333333333 * gamma * curv_weight(p0.x_,p0.y_, p1.x_,p1.y_, x,y,         2.0, bruckstein);
	en += 0.33333333333333333333 * gamma * curv_weight(p1.x_,p1.y_, x,y,         p2.x_,p2.y_, 2.0, bruckstein);
	en += 0.33333333333333333333 * gamma * curv_weight(x,y,         p2.x_,p2.y_, p3.x_,p3.y_, 2.0, bruckstein);

	//en = 0;
	en += 0.5 * integrator.fg_energy_line(p1.x_,p1.y_, x,y);
	en += 0.5 * integrator.fg_energy_line(x,y, p2.x_,p2.y_);

	return en;
}

void SegmentationCurve::draw(std::string svgfile, const Math2D::Matrix<float>& data_term)
{
	size_t xDim = data_term.xDim();
	size_t yDim = data_term.yDim();

	ofstream of(svgfile.c_str());
    init_svg_file(of,uint(xDim),uint(yDim));

	float max_data = data_term.max();
	float min_data = data_term.min();

	//Draw data term background
	for (size_t i=0;i<xDim;++i) {
	for (size_t j=0;j<yDim;++j) {
		vector<pair<double, double> > points;
		points.push_back( make_pair(i,j) );
		points.push_back( make_pair(i+1,j) );
		points.push_back( make_pair(i+1,j+1) );
		points.push_back( make_pair(i,j+1) );
		stringstream sout;
		int c = 255.0*(data_term(i,j)-min_data) / (max_data - min_data);
		sout << "fill:#" << hex << setw(2) << c << c << c
		     << ";stroke:none;";
		svg_draw_polygon(of,sout.str(),points);
	}}


	string line_style = "stroke-width:0.2;stroke:red;";

	//Draw curve
	for (uint i=0; i < original_coord.size(); i++) {
		Mesh2DPoint& p1 = original_coord[i];
		Mesh2DPoint& p2 = original_coord[ (i+1) % original_coord.size() ];
		if (!(p1 == p2)) {
			//Edge still used
			of << "<line style=\"" << line_style << "\"  x1=\"" << p1.x_ << "\" y1=\"" << p1.y_ 
				<< "\" x2=\"" << p2.x_ << "\" y2=\"" << p2.y_ << "\" />" << std::endl;
		}
	}

	line_style = "stroke-width:0.2;stroke:blue;";

	//Draw curve
	for (uint i=0; i < coord.size(); i++) {
		Mesh2DPoint& p1 = coord[i];
		Mesh2DPoint& p2 = coord[ (i+1) % coord.size() ];
		if (!(p1 == p2)) {
			//Edge still used
			of << "<line style=\"" << line_style << "\"  x1=\"" << p1.x_ << "\" y1=\"" << p1.y_ 
				<< "\" x2=\"" << p2.x_ << "\" y2=\"" << p2.y_ << "\" />" << std::endl;
		}
	}
	finish_svg_file(of);
}












