include common/Makefile.common

COMMONDIR = common
CLPLIBDIR = ~/lib/
INCLUDE += -I common/ -I thirdparty/ -I . -I Clp -I ClpInc -I CoinUtils -I CoinUtilsInc -I OsiClp -I Osi  -DHAVE_CONFIG_H -I QPBO-v1.3.src -I HOCR

# if you have the any of the solvers Gurobi, Cplex or Xpress, please add -DHAS_GUROBI etc. to the INCLUDE options

EXTERNALS = QPBO-v1.3.src/opt/*.o thirdparty/gpc.o $(CLPLIBDIR)/libOsiClp.so $(CLPLIBDIR)/libOsi.so $(CLPLIBDIR)/libClp.so $(CLPLIBDIR)/libCoinUtils.so

all: bin/lpseg.debug.L64 bin/lpinpaint.debug.L64 bin/interactiveseg.debug.L64 

bin/lpseg.debug.L64 : lpseg.cc $(DEBUGDIR)/gpcpetter.o $(DEBUGDIR)/lp_segmentation.o $(DEBUGDIR)/extended_lp_segmentation.o $(DEBUGDIR)/lp_segmenter.o $(DEBUGDIR)/segmentation_common.o $(DEBUGDIR)/qpbo_segmentation.o $(DEBUGDIR)/curvature.o $(DEBUGDIR)/mesh2D.o $(DEBUGDIR)/adaptive_mesh.o $(DEBUGDIR)/svg.o $(DEBUGDIR)/draw_segmentation.o $(DEBUGDIR)/conversion.o $(DEBUGDIR)/curvature.o thirdparty/gpc.o common/lib/commonlib.debug 
	$(LINKER) $(DEBUGFLAGS) $(INCLUDE) lpseg.cc $(DEBUGDIR)/gpcpetter.o $(DEBUGDIR)/qpbo_segmentation.o $(DEBUGDIR)/lp_segmentation.o $(DEBUGDIR)/adaptive_mesh.o $(DEBUGDIR)/extended_lp_segmentation.o $(DEBUGDIR)/lp_segmenter.o $(DEBUGDIR)/segmentation_common.o $(DEBUGDIR)/mesh2D.o $(DEBUGDIR)/svg.o $(DEBUGDIR)/draw_segmentation.o $(DEBUGDIR)/conversion.o $(DEBUGDIR)/curvature.o common/lib/commonlib.debug $(EXTERNALS) -lblas -llapack -o $@

bin/interactiveseg.debug.L64 : interactive_seg.cc $(DEBUGDIR)/gpcpetter.o $(DEBUGDIR)/lp_segmentation.o $(DEBUGDIR)/segmentation_common.o $(DEBUGDIR)/curvature.o $(DEBUGDIR)/mesh2D.o $(DEBUGDIR)/svg.o $(DEBUGDIR)/draw_segmentation.o $(DEBUGDIR)/conversion.o $(DEBUGDIR)/curvature.o thirdparty/gpc.o common/lib/commonlib.debug $(DEBUGDIR)/extended_lp_segmentation.o $(DEBUGDIR)/lp_segmenter.o $(DEBUGDIR)/adaptive_mesh.o
	$(LINKER) $(DEBUGFLAGS) $(INCLUDE) interactive_seg.cc $(DEBUGDIR)/gpcpetter.o $(DEBUGDIR)/lp_segmentation.o $(DEBUGDIR)/adaptive_mesh.o $(DEBUGDIR)/segmentation_common.o $(DEBUGDIR)/mesh2D.o $(DEBUGDIR)/svg.o $(DEBUGDIR)/draw_segmentation.o $(DEBUGDIR)/extended_lp_segmentation.o $(DEBUGDIR)/lp_segmenter.o $(DEBUGDIR)/conversion.o $(DEBUGDIR)/curvature.o common/lib/commonlib.debug $(EXTERNALS) -lblas -llapack -o $@

bin/lpinpaint.debug.L64 : lpinpaint.cc $(DEBUGDIR)/gpcpetter.o $(DEBUGDIR)/lp_inpainting.o $(DEBUGDIR)/segmentation_common.o $(DEBUGDIR)/label_components.o $(DEBUGDIR)/color_conversion.o $(DEBUGDIR)/curvature.o $(DEBUGDIR)/mesh2D.o $(DEBUGDIR)/svg.o $(DEBUGDIR)/draw_segmentation.o $(DEBUGDIR)/conversion.o $(DEBUGDIR)/curvature.o thirdparty/gpc.o common/lib/commonlib.debug
	$(LINKER) $(DEBUGFLAGS) $(INCLUDE) lpinpaint.cc $(DEBUGDIR)/gpcpetter.o $(DEBUGDIR)/label_components.o $(DEBUGDIR)/color_conversion.o $(DEBUGDIR)/lp_inpainting.o $(DEBUGDIR)/segmentation_common.o $(DEBUGDIR)/mesh2D.o $(DEBUGDIR)/svg.o $(DEBUGDIR)/draw_segmentation.o $(DEBUGDIR)/conversion.o $(DEBUGDIR)/curvature.o common/lib/commonlib.debug $(EXTERNALS) -lblas -llapack -o $@



include common/Makefile.finish
