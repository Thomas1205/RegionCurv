/*** written by Thomas Schoenemann as a private person without employment, September 2009 ***/

#ifndef SVG_HH
#define SVG_HH

#include <fstream>
#include <vector>
#include "makros.hh"

void init_svg_file(std::ofstream& of, uint max_x, uint max_y);

void finish_svg_file(std::ofstream& of);

void svg_draw_line(std::ofstream& of, std::string style, uint x1, uint y1, uint x2, uint y2);

void svg_draw_rect(std::ofstream& of, std::string style, uint min_x, uint min_y, uint max_x, uint max_y);

void svg_draw_polygon(std::ofstream& of, std::string style, std::vector<std::pair<double,double> >& points);

#endif
