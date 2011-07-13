/*** written by Thomas Schoenemann as a private person without employment, September 2009 ***/

#ifndef SPARSE_MATRIX_DESCRIPTION_HH
#define SPARSE_MATRIX_DESCRIPTION_HH

#include "makros.hh"
#include "vector.hh"
#include <map>

template<typename T>
class SparseMatrixDescription {
public:

  SparseMatrixDescription(uint nEntries, uint nRows = MAX_UINT, uint nColumns = MAX_UINT);

  ~SparseMatrixDescription();

  void sort_by_column(Math1D::Vector<uint>& column_start, bool sort_each_column = false);

  void sort_by_row(Math1D::Vector<uint>& row_start, bool sort_each_row = false);

  const uint* row_indices() const;
  
  const uint* col_indices() const;

  T* value();
  
  uint* row_indices();
  
  uint* col_indices();

  const T* value() const;

  uint nEntries() const;

  uint nReservedEntries() const;

  void add_entry(uint row, uint col, T value);

  //all row- and column indices are inrementeed by one. 
  // this is useful when using GLPK as then indices must start at 1, not at 0 (as e.g. in CLP)
  void increment_indices();

  void increase_nColumns(uint increment);
 
  void reset();

  void reset(uint new_size, int new_nRows = -1, int new_nCols = -1);

protected:

  Math1D::NamedVector<uint> row_idx_;
  Math1D::NamedVector<uint> col_idx_;
  Math1D::NamedVector<T> value_;

  uint nRows_;
  uint nColumns_;
  uint nReservedEntries_;
  uint nActualEntries_;
};



/********************************************* implementation **************************************/

template<typename T>
SparseMatrixDescription<T>::SparseMatrixDescription(uint nEntries, uint nRows, uint nColumns) :
  row_idx_(nEntries,MAKENAME(row_idx_)), col_idx_(nEntries,MAKENAME(col_idx_)),
  value_(nEntries,MAKENAME(value_)),
  nReservedEntries_(nEntries), nActualEntries_(0)  {

  nRows_ = nRows;
  nColumns_ = nColumns;
}

template<typename T>
SparseMatrixDescription<T>::~SparseMatrixDescription() {
}

template<typename T>
void SparseMatrixDescription<T>::increase_nColumns(uint increment) {

  nColumns_ += increment;
}

template<typename T>
void SparseMatrixDescription<T>::sort_by_column(Math1D::Vector<uint>& column_start, bool sort_each_column) {

  int maxColumn = -1;

  std::map<uint,uint> column_count;

  for (uint k=0; k < nActualEntries_; k++) {

    const int col = col_idx_[k];
    if (col > maxColumn)
      maxColumn = col;

    column_count[col]++;
  }

  column_start.resize(maxColumn+2);
  uint s=0;
  for (int k=0; k <= maxColumn; k++) {
    column_start[k] = s;
    s += column_count[k];
  }
  column_start[maxColumn+1] = s;
  
  Math1D::Vector<uint> colcount(maxColumn+1,0);

  Math1D::Vector<uint> temp_row_idx(nReservedEntries_);
  Math1D::Vector<uint> temp_col_idx(nReservedEntries_);
  Math1D::Vector<T> temp_value(nReservedEntries_);

  for (uint k=0; k < nActualEntries_; k++) {

    const uint col = col_idx_[k];
    const uint row = row_idx_[k];
    const T val = value_[k];

    uint idx = column_start[col] + colcount[col];
    colcount[col]++;
    
    temp_row_idx[idx] = row;
    temp_col_idx[idx] = col;
    temp_value[idx] = val;
  }

  if (sort_each_column) {

    for (int col = 0; col < maxColumn; col++) {

      uint size = column_start[col+1] - column_start[col];

      if (size >= 2) {

	//bubble sort
	for (uint k1=1; k1 < size; k1++) {
	  for (uint k2=column_start[col]; k2 < column_start[col+1]-k1; k2++) {
	    
	    if (temp_row_idx[k2+1] < temp_row_idx[k2]) {
	      std::swap(temp_row_idx[k2+1],temp_row_idx[k2]);
	      std::swap(temp_value[k2+1],temp_value[k2]);
	    }
	  }
	}
      }
    }
  }
  
  row_idx_ = temp_row_idx;
  col_idx_ = temp_col_idx;
  value_ = temp_value;
}

template<typename T>
void SparseMatrixDescription<T>::sort_by_row(Math1D::Vector<uint>& row_start, bool sort_each_row) {
  
  std::map<uint,uint> row_count;

  int maxRow = -1;
  for (uint k=0; k < nActualEntries_; k++) {
    const int row = row_idx_[k];
  
    if (maxRow < row)
      maxRow = row;

    row_count[row]++;
  }
  
  row_start.resize(maxRow+2);
  uint s=0;
  for (int k=0; k <= maxRow; k++) {
    row_start[k] = s;
    s += row_count[k];
  }
  row_start[maxRow+1] = s;

  Math1D::Vector<uint> rowcount(maxRow+1,0);

  Math1D::Vector<uint> temp_row_idx(nReservedEntries_);
  Math1D::Vector<uint> temp_col_idx(nReservedEntries_);
  Math1D::Vector<T> temp_value(nReservedEntries_);

  for (uint k=0; k < nActualEntries_; k++) {

    const uint col = col_idx_[k];
    const uint row = row_idx_[k];
    const T val = value_[k];

    uint idx = row_start[row] + rowcount[row];
    rowcount[row]++;
    
    temp_row_idx[idx] = row;
    temp_col_idx[idx] = col;
    temp_value[idx] = val;
  }

  if (sort_each_row) {

    for (int row = 0; row < maxRow; row++) {

      uint size = row_start[row+1] - row_start[row];

      if (size >= 2) {

	for (uint k1=1; k1 < size; k1++) {
	  for (uint k2=row_start[row]; k2 < row_start[row+1]-k1; k2++) {

	    if (temp_col_idx[k2+1] < temp_col_idx[k2]) {
	    
	      std::swap(temp_col_idx[k2+1],temp_col_idx[k2]);
	      std::swap(temp_value[k2+1],temp_value[k2]);
	    }
	  }
	}	
      }
    }
  }
  
  row_idx_ = temp_row_idx;
  col_idx_ = temp_col_idx;
  value_ = temp_value;

}

template<typename T>
void SparseMatrixDescription<T>::reset() {
  nActualEntries_ = 0;
}


template<typename T>
void SparseMatrixDescription<T>::reset(uint new_size, int new_nRows, int new_nCols) {
  nActualEntries_ = 0;
  nReservedEntries_ = new_size;
  row_idx_.resize_dirty(new_size);
  col_idx_.resize_dirty(new_size);
  value_.resize_dirty(new_size);

  if (new_nCols >= 0) {
    nColumns_ = new_nCols;
  }
  if (new_nRows >= 0) {
    nRows_ = new_nRows;
  }
}

template<typename T>
const uint* SparseMatrixDescription<T>::row_indices() const {
  return row_idx_.direct_access();
}

template<typename T>
const uint* SparseMatrixDescription<T>::col_indices() const {
  return col_idx_.direct_access();
}

template<typename T>
const T* SparseMatrixDescription<T>::value() const {
  return value_.direct_access();
}

template<typename T>
uint* SparseMatrixDescription<T>::row_indices() {
  return row_idx_.direct_access();
}

template<typename T>
uint* SparseMatrixDescription<T>::col_indices() {
  return col_idx_.direct_access();
}

template<typename T>
T* SparseMatrixDescription<T>::value() {
  return value_.direct_access();
}

template<typename T>
void SparseMatrixDescription<T>::increment_indices() {

  for (uint i=0; i < nActualEntries_; i++)
    row_idx_[i]++;
  for (uint i=0; i < nActualEntries_; i++)
    col_idx_[i]++;
}

template<typename T>
uint SparseMatrixDescription<T>::nEntries() const {
  return nActualEntries_;
}

template<typename T>
uint SparseMatrixDescription<T>::nReservedEntries() const {
  return nReservedEntries_;
}

template<typename T>
void SparseMatrixDescription<T>::add_entry(uint row, uint col, T value) {

  if (nActualEntries_ >= nReservedEntries_) {
    INTERNAL_ERROR << " number of reserved entries exceeded. Exiting..." << std::endl;
    exit(1);
  }

  if (row >= nRows_) {
    INTERNAL_ERROR << " maximum row exceeded (" << row << ">=" << nRows_ << "). Exiting." << std::endl;
    exit(1);
  }

  if (col >= nColumns_) {
    INTERNAL_ERROR << col << " exceeds the maximum column. Exiting." << std::endl;
    exit(1);
  }

  row_idx_[nActualEntries_] = row;
  col_idx_[nActualEntries_] = col;
  value_[nActualEntries_]   = value;

  nActualEntries_++;
}




#endif
