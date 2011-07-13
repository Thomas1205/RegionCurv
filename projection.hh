/**** written by Thomas Schoenemann as an employee of Lund University, July 2010 ****/

#ifndef PROJECTION_HH
#define PROJECTION_HH

template <typename T>
inline void projection_on_simplex(T* data, const uint nData) {

  /**** reproject on the simplices [Michelot 1986]****/

  assert(nData > 0);
  
  uint nNonZeros = nData;

  //TEST
  uint iter = 0;
  uint start = 0;
  uint end = nData-1;
  //END_TEST

  while (true) {

    //TEST
    iter++;
      
    if ((iter % 5) == 0) {

      while (data[start] == 0.0)
	start++;
      while(data[end] == 0.0)
	end--;
    }
    //END_TEST

    //a) project onto the plane
    T mean_dev = - 1.0;

    //for (uint k=0; k < nData; k++)
    for (uint k=start; k <= end; k++) {
      mean_dev += data[k];
      assert(fabs(data[k] < 1e75));
    }
    
    mean_dev /= nNonZeros;
    assert(!isnan(mean_dev));
      
    //b) subtract mean
    bool all_pos = true;

    const bool first_iter = (nNonZeros == nData);

    //for (uint k=0; k < nData; k++) {
    for (uint k=start; k <= end; k++) {


      T temp = data[k];

      if (first_iter || temp != 0.0) {
	temp -= mean_dev;

	if (temp < 1e-12) {
	  all_pos = false;
	  temp = 0.0;
	  nNonZeros--;
	}
	data[k] = temp;
      }
    }

    if (all_pos)
      break;
  }
}

inline void projection_on_simplex_with_slack(double* data, double& slack, uint nData) {

  uint nNonZeros = nData + 1;

  while (nNonZeros > 0) {
      
    //a) project onto the plane
    double mean_dev = - 1.0 + slack;
    for (uint k=0; k < nData; k++) {
      mean_dev += data[k];
      assert(fabs(data[k] < 1e75));
    }
    
    mean_dev /= nNonZeros;
      
    //b) subtract mean
    bool all_pos = true;
    
    if (nNonZeros == (nData+1) || slack != 0.0) {
      slack -= mean_dev;

      if (slack < 0.0)
	all_pos = false;
    }

    for (uint k=0; k < nData; k++) {
      
      if (nNonZeros == (nData+1) || data[k] != 0.0) {
	data[k] -= mean_dev;
	
	if (data[k] < 0.0)
	  all_pos = false;
      }
    }
    
    if (all_pos)
      break;
    
    //c) fix negatives to 0
    nNonZeros = nData+1;
    if (slack < 1e-8) {
      slack = 0.0;
      nNonZeros--;
    }

    for (uint k=0; k < nData; k++) {
	
      if (data[k] < 1e-8) {
	data[k] = 0.0;
	nNonZeros--;
      }
    }
  }
}

#endif
