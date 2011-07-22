/*** written by Thomas Schoenemann as an employee of Lund University, Sweden, September 2010 ***/


#include "grid.hh"

GridLinePair::GridLinePair(uint line1, uint line2) :
  line1_(line1), line2_(line2) {}

GridLine::GridLine(uint x1, uint y1, uint x2, uint y2) :
  x1_(x1), y1_(y1), x2_(x2), y2_(y2) {}

double GridLine::length() const {

  double dx = ((int) x1_) - ((int) x2_);
  double dy = ((int) y1_) - ((int) y2_);

  return sqrt(dx*dx+dy*dy);
}

Grid::Grid(uint xDim, uint yDim) :
  grid_node_(xDim,yDim) {}


uint Grid::nLines() const {
  return uint( grid_line_.size() );
}

uint Grid::nLinePairs() const {
  return uint( grid_line_pair_.size() );
}

void Grid::add_line(uint x1, uint y1, uint x2, uint y2) {

  uint line_num = grid_line_.size();
  grid_line_.push_back(GridLine(x1,y1,x2,y2));
  
  grid_node_(x1,y1).outgoing_lines_.push_back(line_num);
  grid_node_(x2,y2).incoming_lines_.push_back(line_num);
}

uint Grid::find_line(uint x1, uint y1, uint x2, uint y2) const {


  const std::vector<uint>& outgoing_lines = grid_node_(x1,y1).outgoing_lines_;

  for (uint k=0; k < outgoing_lines.size(); k++) {
    const uint idx = outgoing_lines[k];
    if (grid_line_[idx].x2_ == x2 && grid_line_[idx].y2_ == y2) 
      return idx;
  }

  return MAX_UINT;
}

const GridNode& Grid::get_node(uint x, uint y) const {
  return grid_node_(x,y);
}

const GridLine& Grid::get_line(uint k) const {
  return grid_line_[k];
}

const GridLinePair& Grid::get_line_pair(uint k) const {
  return grid_line_pair_[k];
}

void Grid::list_incoming_line_pairs(uint x, uint y, std::vector<uint>& indices) const {

  indices.clear();

  const std::vector<uint>& in_lines = grid_node_(x,y).incoming_lines_;

  for (uint k=0; k < in_lines.size(); k++) {

    uint l1 = in_lines[k];

    const std::vector<uint>& cur_pairs = grid_line_[l1].ending_line_pairs_ ;

    for (uint k2 = 0; k2 < cur_pairs.size(); k2++) {
      indices.push_back(cur_pairs[k2]);
    }
  }
}

void Grid::list_outgoing_line_pairs(uint x, uint y, std::vector<uint>& indices) const {

  indices.clear();

  const std::vector<uint>& out_lines = grid_node_(x,y).outgoing_lines_;

  for (uint k=0; k < out_lines.size(); k++) {

    uint l1 = out_lines[k];

    const std::vector<uint>& cur_pairs = grid_line_[l1].starting_line_pairs_ ;
    
    for (uint k2 = 0; k2 < cur_pairs.size(); k2++) {
      indices.push_back(cur_pairs[k2]);
    }    
  }  

}



void Grid::generate_grid_line_pairs() {

  for (uint y=0; y < grid_node_.yDim(); y++) {
    for (uint x=0; x < grid_node_.xDim(); x++) {

      const std::vector<uint>& outgoing_lines = grid_node_(x,y).outgoing_lines_;

      for (uint k=0; k < outgoing_lines.size(); k++) {
	
	const uint l1 = outgoing_lines[k];

	const uint hx = grid_line_[l1].x2_;
	const uint hy = grid_line_[l1].y2_;

	//careful: do not list reverse line pairs
	const std::vector<uint>& head_outgoing_lines = grid_node_(hx,hy).outgoing_lines_;

	for (uint k2 = 0; k2 < head_outgoing_lines.size(); k2++) {

	  const uint l2 =head_outgoing_lines[k2];

	  //careful: do not list reverse line pairs
	  if (grid_line_[l2].x2_ != x || grid_line_[l2].y2_ != y) {

	    uint pair_num = uint( grid_line_pair_.size() );
	    grid_line_pair_.push_back(GridLinePair(l1,l2));
	    
	    grid_line_[l1].starting_line_pairs_.push_back(pair_num);
	    grid_line_[l2].ending_line_pairs_.push_back(pair_num);
	  }
	}

      }
    }
  }

}

