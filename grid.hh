/*** written by Thomas Schoenemann as an employee of Lund University, Sweden, September 2010 ***/

#ifndef GRID_HH
#define GRID_HH


#include "vector.hh"
#include "matrix.hh"
#include <vector>

class GridNode {
public:

  std::vector<uint> incoming_lines_;
  std::vector<uint> outgoing_lines_;
};

class GridLine { //note: lines are directed
public:
  
  //note: the line starts at (x1,y1), then proceeds to (x2,y2)
  GridLine(uint x1, uint y1, uint x2, uint y2);

  uint x1_;
  uint y1_;

  uint x2_;
  uint y2_;

  double length() const;

  std::vector<uint> starting_line_pairs_;
  std::vector<uint> ending_line_pairs_;
};

class GridLinePair {
public:

  GridLinePair(uint line1, uint line2);

  uint line1_;
  uint line2_;
};

class Grid {
public:

  Grid(uint xDim, uint yDim);

  void add_line(uint x1, uint y1, uint x2, uint y2);

  uint find_line(uint x1, uint y1, uint x2, uint y2) const;

  uint nLines() const;

  uint nLinePairs() const;

  void list_incoming_line_pairs(uint x, uint y, std::vector<uint>& indices) const;

  void list_outgoing_line_pairs(uint x, uint y, std::vector<uint>& indices) const;

  const GridNode& get_node(uint x, uint y) const;

  const GridLine& get_line(uint k) const;

  const GridLinePair& get_line_pair(uint k) const;

  void generate_grid_line_pairs();

protected:

  Storage2D<GridNode> grid_node_;
  std::vector<GridLine> grid_line_;
  std::vector<GridLinePair> grid_line_pair_;
};



#endif
