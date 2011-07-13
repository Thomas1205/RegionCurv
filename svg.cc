/*** written by Thomas Schoenemann as a private person without employment, September 2009 ***/

#include "svg.hh"

void init_svg_file(std::ofstream& of, uint max_x, uint max_y) {

  assert(of.is_open());

  of << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>" << std::endl;
  of << "<svg width=\"" << max_x << "\" height=\"" << max_y << "\" id=\"svg2\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">" << std::endl; 
}

void finish_svg_file(std::ofstream& of) {

  assert(of.is_open());

  of << "</svg>" << std::endl;
  of.close();
}

void svg_draw_line(std::ofstream& of, std::string style, uint x1, uint y1, uint x2, uint y2) {

  of << "<line style=\"" << style << "\"  x1=\"" << x1 << "\" y1=\"" << y1 
     << "\" x2=\"" << x2 << "\" y2=\"" << y2 << "\" />" << std::endl;
}

void svg_draw_rect(std::ofstream& of, std::string style, uint min_x, uint min_y, uint max_x, uint max_y) {

  of << "<rect style=\"" << style << "\" x=\"" << min_x << "\" y=\"" << min_y 
     << "\" width=\"" << (max_x-min_x) << "\" height=\"" << (max_y-min_y) << "\" />" 
     << std::endl;
}

void svg_draw_polygon(std::ofstream& of, std::string style,
		      std::vector<std::pair<double,double> >& points) {

  of << "<polygon style=\"" << style << "\" points=\"";
  for (std::vector<std::pair<double,double> >::iterator it = points.begin(); it != points.end(); it++)
    of << it->first << "," << it->second << " ";
  of << "\" />" << std::endl;
}


