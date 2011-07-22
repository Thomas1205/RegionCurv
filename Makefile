include common/Makefile.common


# Change these to point to your Coin Cbc installation
COININCLUDEDIR = ~/Programming/coin-Cbc-win/build/include/
COINLIBDIR = ~/Programming/coin-Cbc-win/build/lib/

##############################################################

COMMONDIR = common
INCLUDE = -I common/ -I thirdparty/ -I . -I $(COININCLUDEDIR) -DHAVE_CONFIG_H -I QPBO-v1.3.src -I HOCR

# if you have the any of the solvers Gurobi, Cplex or Xpress, please add -DHAS_GUROBI etc. to the INCLUDE options

QPBO = QPBO-v1.3.src/QPBO.o QPBO-v1.3.src/QPBO_extra.o QPBO-v1.3.src/QPBO_maxflow.o QPBO-v1.3.src/QPBO_postprocessing.o 
GPC = thirdparty/gpc.o
EXTERNALS = -lClp -lCoinUtils -lOsiClp -lOsi -lCbc -lCgl

all: bin/lpseg.debug.L64 bin/lpinpaint.debug.L64 bin/interactiveseg.debug.L64 

bin/lpseg.debug.L64 : lpseg.cc $(DEBUGDIR)/gpcpetter.o $(DEBUGDIR)/lp_segmentation.o $(DEBUGDIR)/extended_lp_segmentation.o $(DEBUGDIR)/lp_segmenter.o $(DEBUGDIR)/segmentation_common.o $(DEBUGDIR)/qpbo_segmentation.o $(DEBUGDIR)/curvature.o $(DEBUGDIR)/mesh2D.o $(DEBUGDIR)/grid.o $(DEBUGDIR)/adaptive_mesh.o $(DEBUGDIR)/svg.o $(DEBUGDIR)/draw_segmentation.o $(DEBUGDIR)/conversion.o $(DEBUGDIR)/curvature.o $(QPBO) $(GPC) common/lib/commonlib.debug 
	$(LINKER) $(DEBUGFLAGS) -L $(COINLIBDIR) $(INCLUDE) lpseg.cc $(DEBUGDIR)/gpcpetter.o $(DEBUGDIR)/qpbo_segmentation.o $(DEBUGDIR)/lp_segmentation.o $(DEBUGDIR)/adaptive_mesh.o $(DEBUGDIR)/extended_lp_segmentation.o $(DEBUGDIR)/lp_segmenter.o $(DEBUGDIR)/segmentation_common.o $(DEBUGDIR)/mesh2D.o $(DEBUGDIR)/grid.o $(DEBUGDIR)/svg.o $(DEBUGDIR)/draw_segmentation.o $(DEBUGDIR)/conversion.o $(DEBUGDIR)/curvature.o common/lib/commonlib.debug $(QPBO) $(GPC) $(EXTERNALS) -lblas -llapack -o $@

bin/interactiveseg.debug.L64 : interactive_seg.cc $(DEBUGDIR)/gpcpetter.o $(DEBUGDIR)/lp_segmentation.o $(DEBUGDIR)/segmentation_common.o $(DEBUGDIR)/curvature.o $(DEBUGDIR)/mesh2D.o $(DEBUGDIR)/grid.o $(DEBUGDIR)/svg.o $(DEBUGDIR)/draw_segmentation.o $(DEBUGDIR)/conversion.o $(DEBUGDIR)/curvature.o $(QPBO) $(GPC) common/lib/commonlib.debug $(DEBUGDIR)/extended_lp_segmentation.o $(DEBUGDIR)/lp_segmenter.o $(DEBUGDIR)/adaptive_mesh.o
	$(LINKER) $(DEBUGFLAGS) -L $(COINLIBDIR) $(INCLUDE) interactive_seg.cc $(DEBUGDIR)/gpcpetter.o $(DEBUGDIR)/lp_segmentation.o $(DEBUGDIR)/adaptive_mesh.o $(DEBUGDIR)/grid.o $(DEBUGDIR)/segmentation_common.o $(DEBUGDIR)/mesh2D.o $(DEBUGDIR)/svg.o $(DEBUGDIR)/draw_segmentation.o $(DEBUGDIR)/extended_lp_segmentation.o $(DEBUGDIR)/lp_segmenter.o $(DEBUGDIR)/conversion.o $(DEBUGDIR)/curvature.o common/lib/commonlib.debug $(QPBO) $(GPC) $(EXTERNALS) -lblas -llapack -o $@

bin/lpinpaint.debug.L64 : lpinpaint.cc $(DEBUGDIR)/gpcpetter.o $(DEBUGDIR)/lp_inpainting.o $(DEBUGDIR)/segmentation_common.o $(DEBUGDIR)/label_components.o $(DEBUGDIR)/color_conversion.o $(DEBUGDIR)/curvature.o $(DEBUGDIR)/mesh2D.o $(DEBUGDIR)/grid.o $(DEBUGDIR)/svg.o $(DEBUGDIR)/draw_segmentation.o $(DEBUGDIR)/conversion.o $(DEBUGDIR)/curvature.o $(QPBO) $(GPC) common/lib/commonlib.debug
	$(LINKER) -L $(COINLIBDIR) $(DEBUGFLAGS) $(INCLUDE) lpinpaint.cc $(DEBUGDIR)/gpcpetter.o $(DEBUGDIR)/label_components.o $(DEBUGDIR)/color_conversion.o $(DEBUGDIR)/lp_inpainting.o $(DEBUGDIR)/segmentation_common.o $(DEBUGDIR)/mesh2D.o $(DEBUGDIR)/grid.o $(DEBUGDIR)/svg.o $(DEBUGDIR)/draw_segmentation.o $(DEBUGDIR)/conversion.o $(DEBUGDIR)/curvature.o common/lib/commonlib.debug $(QPBO) $(GPC) $(EXTERNALS) -lblas -llapack -o $@



include common/Makefile.finish
