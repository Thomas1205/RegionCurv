/******** written by Thomas Schoenemann as an employee of Lund University, Sweden, October 2010 ****/

#ifndef CURV_DENOISING
#define CURV_DENOISING

#include "matrix.hh"
#include "segmentation_common.hh"

//NOTE: currently only the first channel is denoised
void curv_denoise(const Math3D::Tensor<float>& image, const LPSegOptions& options, 
		  Math3D::Tensor<double>& denoised_image, uint nBins=1);


#endif
