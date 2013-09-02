include common/Makefile.common


# Change these to point to your Coin Cbc installation
COININCLUDEDIR = ~/work/software/linprog/thirdparty/devel/coin-Cbc/include/
COINLIBDIR = ~/work/software/linprog/thirdparty/devel/coin-Cbc/lib64/

##############################################################

COMMONDIR = common
INCLUDE = -I common/ -I msg_passing/ -I thirdparty/ -I . -I $(COININCLUDEDIR) -DHAVE_CONFIG_H -I QPBO-v1.3.src -I HOCR

# if you have the any of the solvers Gurobi, Cplex or Xpress, please add -DHAS_GUROBI etc. to the INCLUDE options

#if you don't have QPBO, outcomment the next three lines

DEBUGFLAGS += -DHAS_QPBO
OPTFLAGS += -DHAS_QPBO
QPBO = QPBO-v1.3.src/QPBO.o QPBO-v1.3.src/QPBO_extra.o QPBO-v1.3.src/QPBO_maxflow.o QPBO-v1.3.src/QPBO_postprocessing.o 


GPC = thirdparty/gpc.o
EXTERNALS = -lClp -lCoinUtils -lOsiClp -lOsi -lCbc -lCgl -lz -lbz2

all: .common bin $(DEBUGDIR) $(OPTDIR) bin/lpseg.opt.L64 bin/lpinpaint.debug.L64 bin/interactiveseg.opt.L64 bin/curvdenoise.debug.L64 common/lib/commonlib.debug common/lib/commonlib.opt

.common: 
	cd common; make; cd -

common/lib/commonlib.opt : 
	cd common; make; cd -

common/lib/commonlib.debug :
	cd common; make; cd -

bin/lpseg.opt.L64 : lpseg.cc $(OPTDIR)/gpcpetter.o $(OPTDIR)/lp_segmentation.o $(OPTDIR)/extended_lp_segmentation.o $(OPTDIR)/lp_segmenter.o $(OPTDIR)/segmentation_common.o $(OPTDIR)/qpbo_segmentation.o $(OPTDIR)/curvature.o $(OPTDIR)/mesh2D.o $(OPTDIR)/grid.o $(OPTDIR)/adaptive_mesh.o $(OPTDIR)/svg.o $(OPTDIR)/draw_segmentation.o $(OPTDIR)/conversion.o $(OPTDIR)/curvature.o $(OPTDIR)/ImageIntegrator.o $(QPBO) $(GPC) common/lib/commonlib.debug msg_passing/$(OPTDIR)/factorChainDualDecomp.o msg_passing/$(OPTDIR)/factorDualOpt.o msg_passing/$(OPTDIR)/factorMPBP.o msg_passing/$(OPTDIR)/factorTRWS.o msg_passing/$(OPTDIR)/separatorChainDualDecomp.o msg_passing/$(OPTDIR)/separatorDualOpt.o msg_passing/$(OPTDIR)/sepTRWS.o
	$(LINKER) $(OPTFLAGS) -L $(COINLIBDIR) $(INCLUDE) lpseg.cc $(OPTDIR)/gpcpetter.o $(OPTDIR)/qpbo_segmentation.o $(OPTDIR)/lp_segmentation.o $(OPTDIR)/adaptive_mesh.o $(OPTDIR)/extended_lp_segmentation.o $(OPTDIR)/lp_segmenter.o $(OPTDIR)/segmentation_common.o $(OPTDIR)/mesh2D.o $(OPTDIR)/grid.o $(OPTDIR)/svg.o $(OPTDIR)/draw_segmentation.o $(OPTDIR)/conversion.o $(OPTDIR)/curvature.o $(OPTDIR)/ImageIntegrator.o common/lib/commonlib.debug msg_passing/$(OPTDIR)/factorChainDualDecomp.o msg_passing/$(OPTDIR)/factorDualOpt.o msg_passing/$(OPTDIR)/factorMPBP.o msg_passing/$(OPTDIR)/factorTRWS.o msg_passing/$(OPTDIR)/separatorChainDualDecomp.o msg_passing/$(OPTDIR)/separatorDualOpt.o msg_passing/$(OPTDIR)/sepTRWS.o $(QPBO) $(GPC) $(EXTERNALS) -lblas -llapack -o $@

bin/interactiveseg.opt.L64 : interactive_seg.cc $(OPTDIR)/gpcpetter.o $(OPTDIR)/lp_segmentation.o $(OPTDIR)/segmentation_common.o $(OPTDIR)/curvature.o $(OPTDIR)/mesh2D.o $(OPTDIR)/grid.o $(OPTDIR)/svg.o $(OPTDIR)/draw_segmentation.o $(OPTDIR)/conversion.o $(OPTDIR)/curvature.o $(QPBO) $(GPC) common/lib/commonlib.debug $(OPTDIR)/extended_lp_segmentation.o $(OPTDIR)/lp_segmenter.o $(OPTDIR)/adaptive_mesh.o $(OPTDIR)/ImageIntegrator.o msg_passing/$(OPTDIR)/factorChainDualDecomp.o msg_passing/$(OPTDIR)/factorDualOpt.o msg_passing/$(OPTDIR)/factorMPBP.o msg_passing/$(OPTDIR)/factorTRWS.o msg_passing/$(OPTDIR)/separatorChainDualDecomp.o msg_passing/$(OPTDIR)/separatorDualOpt.o msg_passing/$(OPTDIR)/sepTRWS.o
	$(LINKER) $(OPTFLAGS) -L $(COINLIBDIR) $(INCLUDE) interactive_seg.cc $(OPTDIR)/gpcpetter.o $(OPTDIR)/lp_segmentation.o $(OPTDIR)/adaptive_mesh.o $(OPTDIR)/grid.o $(OPTDIR)/segmentation_common.o $(OPTDIR)/mesh2D.o $(OPTDIR)/svg.o $(OPTDIR)/draw_segmentation.o $(OPTDIR)/extended_lp_segmentation.o $(OPTDIR)/lp_segmenter.o $(OPTDIR)/conversion.o $(OPTDIR)/curvature.o $(OPTDIR)/ImageIntegrator.o common/lib/commonlib.debug msg_passing/$(OPTDIR)/factorChainDualDecomp.o msg_passing/$(OPTDIR)/factorDualOpt.o msg_passing/$(OPTDIR)/factorMPBP.o msg_passing/$(OPTDIR)/factorTRWS.o msg_passing/$(OPTDIR)/separatorChainDualDecomp.o msg_passing/$(OPTDIR)/separatorDualOpt.o msg_passing/$(OPTDIR)/sepTRWS.o common/$(OPTDIR)/storage2D.o $(QPBO) $(GPC) $(EXTERNALS) -lblas -llapack -o $@

bin/lpinpaint.debug.L64 : lpinpaint.cc $(DEBUGDIR)/gpcpetter.o $(DEBUGDIR)/lp_inpainting.o $(DEBUGDIR)/segmentation_common.o $(DEBUGDIR)/adaptive_mesh.o $(DEBUGDIR)/label_components.o $(DEBUGDIR)/color_conversion.o $(DEBUGDIR)/curvature.o $(DEBUGDIR)/mesh2D.o $(DEBUGDIR)/grid.o $(DEBUGDIR)/svg.o $(DEBUGDIR)/draw_segmentation.o $(DEBUGDIR)/conversion.o $(DEBUGDIR)/curvature.o $(QPBO) $(GPC) common/lib/commonlib.debug
	$(LINKER) -L $(COINLIBDIR) $(DEBUGFLAGS) $(INCLUDE) lpinpaint.cc $(DEBUGDIR)/gpcpetter.o $(DEBUGDIR)/label_components.o $(DEBUGDIR)/color_conversion.o $(DEBUGDIR)/lp_inpainting.o $(DEBUGDIR)/segmentation_common.o $(DEBUGDIR)/adaptive_mesh.o $(DEBUGDIR)/mesh2D.o $(DEBUGDIR)/grid.o $(DEBUGDIR)/svg.o $(DEBUGDIR)/draw_segmentation.o $(DEBUGDIR)/conversion.o $(DEBUGDIR)/curvature.o common/lib/commonlib.debug $(QPBO) $(GPC) $(EXTERNALS) -lblas -llapack -o $@

bin/curvdenoise.debug.L64 : curvdenoise.cc $(DEBUGDIR)/gpcpetter.o $(DEBUGDIR)/curv_denoising.o $(DEBUGDIR)/segmentation_common.o $(DEBUGDIR)/adaptive_mesh.o $(DEBUGDIR)/label_components.o $(DEBUGDIR)/color_conversion.o $(DEBUGDIR)/curvature.o $(DEBUGDIR)/mesh2D.o $(DEBUGDIR)/grid.o $(DEBUGDIR)/svg.o $(DEBUGDIR)/draw_segmentation.o $(DEBUGDIR)/conversion.o $(DEBUGDIR)/curvature.o $(QPBO) $(GPC) common/lib/commonlib.debug
	$(LINKER) -L $(COINLIBDIR) $(DEBUGFLAGS) $(INCLUDE) curvdenoise.cc $(DEBUGDIR)/gpcpetter.o $(DEBUGDIR)/label_components.o $(DEBUGDIR)/adaptive_mesh.o $(DEBUGDIR)/color_conversion.o $(DEBUGDIR)/curv_denoising.o $(DEBUGDIR)/segmentation_common.o $(DEBUGDIR)/mesh2D.o $(DEBUGDIR)/grid.o $(DEBUGDIR)/svg.o $(DEBUGDIR)/draw_segmentation.o $(DEBUGDIR)/conversion.o $(DEBUGDIR)/curvature.o common/lib/commonlib.debug $(QPBO) $(GPC) $(EXTERNALS) -lblas -llapack -o $@

include common/Makefile.finish
