project(superlu_mt)

if (SUPERLU_MT_SOURCE_FOUND)

set(ALLAUX superlu_timer.c dclock.c sp_ienv.c lsame.c xerbla.c 
	util.c pmemory.c qrnzcnt.c await.c 
	get_perm_c.c mmd.c colamd.c sp_coletree.c 
	pxgstrf_scheduler.c sp_colorder.c 
	pxgstrf_mark_busy_descends.c 
	pxgstrf_pruneL.c pxgstrf_super_bnd_dfs.c 
	pxgstrf_relax_snode.c heap_relax_snode.c 
	pxgstrf_synch.c pxgstrf_finalize.c)

set(SZAUX  slamch.c)
set(DZAUX  dlamch.c)

set(SLUSRC sreadhb.c smatgen.c 
	psgssv.c psgssvx.c sgstrs.c sgsrfs.c 
	sgscon.c slacon.c slangs.c sgsequ.c slaqgs.c spivotgrowth.c 
	psmemory.c psutil.c smyblas2.c 
	psgstrf.c psgstrf_init.c 
	psgstrf_thread.c psgstrf_thread_init.c psgstrf_thread_finalize.c 
	psgstrf_factor_snode.c psgstrf_snode_dfs.c psgstrf_snode_bmod.c 
	psgstrf_panel_dfs.c psgstrf_panel_bmod.c 
	psgstrf_copy_to_ucol.c psgstrf_pivotL.c 
	psgstrf_column_dfs.c psgstrf_column_bmod.c 
	psgstrf_bmod1D.c psgstrf_bmod2D.c 
	psgstrf_bmod1D_mv2.c psgstrf_bmod2D_mv2.c 
	ssp_blas2.c ssp_blas3.c)

set(DLUSRC dreadhb.c dmatgen.c 
	pdgssv.c pdgssvx.c dgstrs.c dgsrfs.c 
	dgscon.c dlacon.c dlangs.c dgsequ.c dlaqgs.c dpivotgrowth.c 
	pdmemory.c pdutil.c dmyblas2.c 
	pdgstrf.c pdgstrf_init.c 
	pdgstrf_thread.c pdgstrf_thread_init.c pdgstrf_thread_finalize.c 
	pdgstrf_factor_snode.c pdgstrf_snode_dfs.c pdgstrf_snode_bmod.c 
	pdgstrf_panel_dfs.c pdgstrf_panel_bmod.c 
	pdgstrf_copy_to_ucol.c pdgstrf_pivotL.c 
	pdgstrf_column_dfs.c pdgstrf_column_bmod.c 
	pdgstrf_bmod1D.c pdgstrf_bmod2D.c 
	pdgstrf_bmod1D_mv2.c pdgstrf_bmod2D_mv2.c 
	dsp_blas2.c dsp_blas3.c) 

set(CLUSRC scomplex.c creadhb.c cmatgen.c scsum1.c icmax1.c 
	pcgssv.c pcgssvx.c cgstrs.c cgsrfs.c 
	cgscon.c clacon.c clangs.c cgsequ.c claqgs.c cpivotgrowth.c 
	pcmemory.c pcutil.c cmyblas2.c 
	pcgstrf.c pcgstrf_init.c 
	pcgstrf_thread.c pcgstrf_thread_init.c pcgstrf_thread_finalize.c 
	pcgstrf_factor_snode.c pcgstrf_snode_dfs.c pcgstrf_snode_bmod.c 
	pcgstrf_panel_dfs.c pcgstrf_panel_bmod.c 
	pcgstrf_copy_to_ucol.c pcgstrf_pivotL.c 
	pcgstrf_column_dfs.c pcgstrf_column_bmod.c 
	pcgstrf_bmod1D.c pcgstrf_bmod2D.c 
	pcgstrf_bmod1D_mv2.c pcgstrf_bmod2D_mv2.c 
	csp_blas2.c csp_blas3.c)

set(ZLUSRC dcomplex.c zreadhb.c zmatgen.c dzsum1.c izmax1.c 
	pzgssv.c pzgssvx.c zgstrs.c zgsrfs.c 
	zgscon.c zlacon.c zlangs.c zgsequ.c zlaqgs.c zpivotgrowth.c 
	pzmemory.c pzutil.c zmyblas2.c 
	pzgstrf.c pzgstrf_init.c 
	pzgstrf_thread.c pzgstrf_thread_init.c pzgstrf_thread_finalize.c 
	pzgstrf_factor_snode.c pzgstrf_snode_dfs.c pzgstrf_snode_bmod.c 
	pzgstrf_panel_dfs.c pzgstrf_panel_bmod.c 
	pzgstrf_copy_to_ucol.c pzgstrf_pivotL.c 
	pzgstrf_column_dfs.c pzgstrf_column_bmod.c 
	pzgstrf_bmod1D.c pzgstrf_bmod2D.c 
	pzgstrf_bmod1D_mv2.c pzgstrf_bmod2D_mv2.c 
	zsp_blas2.c zsp_blas3.c)

add_definitions(-D__OPENMP)
add_library(superlu_mt ${ALLAUX} ${SLUSRC} ${DLUSRC} ${CLUSRC} ${ZLUSRC})

endif (SUPERLU_MT_SOURCE_FOUND)
