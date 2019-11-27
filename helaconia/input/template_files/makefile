
FC= gfortran
objdir=$(PWD)/obj/
moddir=$(PWD)/mod/
srcdir=$(PWD)/src/
decaydir=$(PWD)/decay/
libdir=$(PWD)/lib/
plotdir=$(PWD)/analysis/user/
showerdir=$(PWD)/shower/interface/
lhapdfobj= $(objdir)Structf_LHAPDFs.o 
obj=$(objdir)Helac_Global.o $(objdir)Helac_Func_1.o  $(objdir)Feynman_Helac.o $(objdir)SM_FeynRule_Helac.o $(objdir)Constants.o $(objdir)Structf_PDFs.o $(objdir)CTEQ6PDF.o $(objdir)Helac_histo.o $(objdir)Helac_master.o $(objdir)Helac_pan1.o $(objdir)Helac_pan2.o $(objdir)Helac_wavef.o $(objdir)Helac_ranmar.o $(objdir)Projectors.o $(objdir)ADAPT.o $(objdir)Kinetic_Func.o $(objdir)KT_Clustering.o $(objdir)MC_Funcs.o $(objdir)MC_Helac_GRID.o $(objdir)MC_PARNI_Weight.o $(objdir)MC_RAMBO.o $(objdir)unweight_lhe.o $(objdir)Phegas_Durham.o $(objdir)Phegas_Choice.o $(objdir)Phegas.o $(objdir)Cuts_Module.o $(objdir)Helac_unwei.o $(objdir)SinglePro.o $(objdir)Summation_Pro.o $(objdir)Main_Test.o $(objdir)MC_VEGAS.o $(objdir)Func_PSI.o $(objdir)Colliders_PSI1.o $(objdir)Colliders_PSI2.o $(objdir)setscale.o $(objdir)alfas_functions.o $(lhapdfobj) $(objdir)FO_plot.o $(objdir)plot_user.o $(objdir)QEDPS_interface.o $(objdir)mint-integrator.o $(objdir)DecayInfo.o $(objdir)HOVll.o $(objdir)Decay_interface.o $(objdir)reweight_xsec.o $(objdir)HO_t2bw.o $(objdir)HO_chi2psia.o
mod = $(moddir)/*.mod
exe = $(PWD)/bin/Helac-Onia
linkexe = $(PWD)/Helac-Onia
EXLIBS =  -L/Users/erdissshaw/Works/LHAPDF/lhapdf_installing/lib -lLHAPDF
MAKEMORE =  makefile_fjcore makefile_fjcore_heptoptagger
all :   clean $(exe)
	ln -s $(exe) $(linkexe)

$(exe): makefile_pdf make_qedps make_hbook make_minuit $(MAKEMORE) $(obj)
	@cp *.o $(objdir)
	@cp *.mod $(moddir)
	@rm *.o *.mod
	$(FC) -o $(exe) $(obj) -L$(libdir) -lpdf -lqedps -lhbook -lminuit $(EXLIBS)

$(objdir)Helac_Global.o $(moddir)Helac_Global.mod : $(srcdir)Helac_Global.f90
	$(FC) -c $< $(FFLAGS)
$(objdir)Helac_Func_1.o $(moddir)Helac_Func_1.mod : $(srcdir)Helac_Func_1.f90 $(moddir)Helac_Global.mod $(moddir)alfas_functions.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)Feynman_Helac.o $(moddir)Feynman_Helac.mod: $(srcdir)Feynman_Helac.f90 $(moddir)Helac_Func_1.mod $(moddir)Helac_Global.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)SM_FeynRule_Helac.o $(moddir)SM_FeynRule_Helac.mod: $(srcdir)SM_FeynRule_Helac.f90 $(moddir)Helac_Global.mod $(moddir)Helac_Func_1.mod $(moddir)Structf_PDFs.mod $(moddir)Constants.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)Constants.o $(moddir)Constants.mod: $(srcdir)Constants.f90 $(moddir)Helac_Global.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)Structf_PDFs.o $(moddir)Structf_PDFs.mod:$(srcdir)Structf_PDFs.f90 $(moddir)Helac_Global.mod $(moddir)CTEQ6PDF.mod $(moddir)setscale.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)CTEQ6PDF.o $(moddir)CTEQ6PDF.mod: $(srcdir)CTEQ6PDF.f90
	$(FC) -c $< $(FFLAGS)
$(objdir)Helac_histo.o $(moddir)Helac_histo.mod:$(srcdir)Helac_histo.f90
	$(FC) -c $< $(FFLAGS)
$(objdir)Helac_master.o $(moddir)Helac_master.mod: $(srcdir)Helac_master.f90 $(moddir)Helac_Global.mod $(moddir)Helac_Func_1.mod $(moddir)Helac_pan1.mod $(moddir)Feynman_Helac.mod $(moddir)Projectors.mod $(moddir)Decay_interface.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)Helac_pan1.o $(moddir)Helac_pan1.mod:$(srcdir)Helac_pan1.f90 $(moddir)Helac_Global.mod $(moddir)Helac_Func_1.mod $(moddir)Feynman_Helac.mod $(moddir)Helac_wavef.mod $(moddir)SM_FeynRule_Helac.mod $(moddir)Helac_pan2.mod $(moddir)Helac_ranmar.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)Helac_pan2.o $(moddir)Helac_pan2.mod: $(srcdir)Helac_pan2.f90 $(moddir)Helac_Global.mod $(moddir)Helac_Func_1.mod $(moddir)Helac_wavef.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)Helac_wavef.o $(moddir)Helac_wavef.mod: $(srcdir)Helac_wavef.f90 $(moddir)Helac_Global.mod $(moddir)Helac_Func_1.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)Helac_ranmar.o $(moddir)Helac_ranmar.mod: $(srcdir)Helac_ranmar.f90
	$(FC) -c $< $(FFLAGS)
$(objdir)Projectors.o $(moddir)Projectors.mod: $(srcdir)Projectors.f90 $(moddir)Helac_wavef.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)ADAPT.o $(moddir)ADAPT.mod: $(srcdir)ADAPT.f90 $(moddir)Helac_Global.mod $(moddir)Helac_ranmar.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)Kinetic_Func.o $(moddir)Kinetic_Func.mod: $(srcdir)Kinetic_Func.f90 $(moddir)Helac_Global.mod $(moddir)Helac_Func_1.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)KT_Clustering.o $(moddir)KT_Clustering.mod:  $(srcdir)KT_Clustering.f90 $(moddir)Helac_Global.mod $(moddir)Structf_PDFs.mod $(moddir)Helac_Func_1.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)MC_Funcs.o $(moddir)MC_Funcs.mod:  $(srcdir)MC_Funcs.f90 $(moddir)Helac_Global.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)MC_Helac_GRID.o $(moddir)MC_Helac_GRID.mod:  $(srcdir)MC_Helac_GRID.f90 $(moddir)Helac_ranmar.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)MC_PARNI_Weight.o $(moddir)MC_PARNI_Weight.mod:  $(srcdir)MC_PARNI_Weight.f90
	$(FC) -c $< $(FFLAGS)
$(objdir)MC_RAMBO.o $(moddir)MC_RAMBO.mod: $(srcdir)MC_RAMBO.f90 $(moddir)Helac_ranmar.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)unweight_lhe.o $(moddir)unweight_lhe.mod:  $(srcdir)unweight_lhe.f90
	$(FC) -c $< $(FFLAGS)
$(objdir)Phegas_Durham.o $(moddir)Phegas_Durham.mod:  $(srcdir)Phegas_Durham.f90 $(moddir)Helac_Func_1.mod $(moddir)Kinetic_Func.mod $(moddir)Helac_ranmar.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)Phegas_Choice.o $(moddir)Phegas_Choice.mod: $(srcdir)Phegas_Choice.f90
	$(FC) -c $< $(FFLAGS)
$(objdir)Phegas.o $(moddir)Phegas.mod: $(srcdir)Phegas.f90 $(moddir)Helac_Global.mod $(moddir)Helac_Func_1.mod $(moddir)Helac_pan2.mod $(moddir)Kinetic_Func.mod $(moddir)MC_RAMBO.mod $(moddir)Phegas_Durham.mod $(moddir)Feynman_Helac.mod $(moddir)Phegas_Choice.mod $(moddir)MC_Helac_GRID.mod $(moddir)Cuts_Module.mod $(moddir)Helac_histo.mod $(moddir)ADAPT.mod $(moddir)QEDPS_interface.mod $(moddir)Decay_interface.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)Cuts_Module.o $(moddir)Cuts_Module.mod: $(srcdir)Cuts_Module.f90 $(moddir)Helac_Global.mod $(moddir)Kinetic_Func.mod $(moddir)Constants.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)Helac_unwei.o $(moddir)Helac_unwei.mod: $(srcdir)Helac_unwei.f90 $(moddir)Helac_histo.mod $(moddir)Helac_Global.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)SinglePro.o $(moddir)SinglePro.mod: $(srcdir)SinglePro.f90 $(moddir)Helac_Global.mod $(moddir)Helac_master.mod $(moddir)Helac_histo.mod $(moddir)Constants.mod $(moddir)unweight_lhe.mod $(moddir)Feynman_Helac.mod $(moddir)Helac_Func_1.mod $(moddir)Phegas.mod $(moddir)Cuts_Module.mod $(moddir)SM_FeynRule_Helac.mod $(moddir)Helac_ranmar.mod $(moddir)Kinetic_Func.mod $(moddir)MC_Funcs.mod $(moddir)Structf_PDFs.mod $(moddir)KT_Clustering.mod $(moddir)Helac_unwei.mod $(moddir)MC_Helac_GRID.mod $(moddir)ADAPT.mod $(moddir)Colliders_PSI1.mod $(moddir)Colliders_PSI2.mod $(moddir)FO_plot.mod $(moddir)QEDPS_interface.mod $(moddir)DecayInfo.mod $(moddir)Decay_interface.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)Summation_Pro.o $(moddir)Summation_Pro.mod: $(srcdir)Summation_Pro.f90 $(moddir)Helac_Global.mod $(moddir)Constants.mod $(moddir)SinglePro.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)Main_Test.o: $(srcdir)Main_Test.f90 $(moddir)SinglePro.mod $(moddir)Summation_Pro.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)MC_VEGAS.o $(moddir)MC_VEGAS.mod: $(srcdir)MC_VEGAS.f90
	$(FC) -c $< $(FFLAGS)
$(objdir)Func_PSI.o $(moddir)Func_PSI.mod: $(srcdir)Func_PSI.f90
	$(FC) -c $< $(FFLAGS)
$(objdir)Colliders_PSI1.o $(moddir)Colliders_PSI1.mod: $(srcdir)Colliders_PSI1.f90 $(moddir)MC_VEGAS.mod $(moddir)Helac_Global.mod $(moddir)Constants.mod $(moddir)Func_PSI.mod $(moddir)SM_FeynRule_Helac.mod $(moddir)Structf_PDFs.mod $(moddir)Helac_master.mod $(moddir)Helac_Func_1.mod $(moddir)Cuts_Module.mod $(moddir)FO_plot.mod $(moddir)QEDPS_interface.mod $(moddir)Decay_interface.mod $(moddir)reweight_xsec.mod $(moddir)KT_Clustering.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)Colliders_PSI2.o $(moddir)Colliders_PSI2.mod: $(srcdir)Colliders_PSI2.f90 $(moddir)Colliders_PSI1.mod $(moddir)MC_VEGAS.mod $(moddir)Helac_Global.mod $(moddir)Constants.mod $(moddir)Func_PSI.mod $(moddir)SM_FeynRule_Helac.mod $(moddir)Structf_PDFs.mod $(moddir)Helac_master.mod $(moddir)Helac_Func_1.mod $(moddir)Cuts_Module.mod $(moddir)FO_plot.mod $(moddir)QEDPS_interface.mod $(moddir)Decay_interface.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)alfas_functions.o $(moddir)alfas_functions.mod: $(srcdir)alfas_functions.f90 $(moddir)Helac_Global.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)setscale.o $(moddir)setscale.mod: $(srcdir)setscale.f90
	$(FC) -c $< $(FFLAGS)
$(objdir)Structf_LHAPDFs.o : $(srcdir)Structf_LHAPDFs.f90 $(moddir)Helac_Global.mod $(moddir)setscale.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)FO_plot.o $(moddir)FO_plot.mod : $(srcdir)FO_plot.f90 $(moddir)Helac_Global.mod $(moddir)plot_user.mod $(moddir)MC_VEGAS.mod $(moddir)Helac_Func_1.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)plot_user.o $(moddir)plot_user.mod : $(plotdir)plot_user.f90 $(moddir)Helac_Global.mod $(moddir)Kinetic_Func.mod $(moddir)QEDPS_interface.mod $(moddir)Decay_interface.mod $(moddir)Helac_Func_1.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)QEDPS_interface.o $(moddir)QEDPS_interface.mod : $(showerdir)QEDPS_interface.f90 $(moddir)Helac_Global.mod $(moddir)Helac_Func_1.mod $(moddir)Kinetic_Func.mod $(moddir)setscale.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)mint-integrator.o $(moddir)mint-integrator.mod : $(srcdir)mint-integrator.f90 $(moddir)Helac_ranmar.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)DecayInfo.o $(moddir)DecayInfo.mod : $(decaydir)DecayInfo.f90 $(moddir)Helac_Global.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)HOVll.o $(moddir)HOVll.mod : $(decaydir)HOVll.f90 $(moddir)Helac_ranmar.mod $(moddir)Kinetic_Func.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)HO_t2bw.o $(moddir)HO_t2bw.mod : $(decaydir)HO_t2bw.f90 $(moddir)Helac_ranmar.mod $(moddir)Kinetic_Func.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)Decay_interface.o $(moddir)Decay_interface.mod : $(decaydir)Decay_interface.f90 $(moddir)Helac_Global.mod $(moddir)HOVll.mod $(moddir)Helac_Func_1.mod $(moddir)Kinetic_Func.mod $(moddir)QEDPS_interface.mod $(moddir)Constants.mod $(moddir)Helac_ranmar.mod $(moddir)HO_t2bw.mod $(moddir)HO_chi2psia.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)reweight_xsec.o $(moddir)reweight_xsec.mod : $(srcdir)reweight_xsec.f90 $(moddir)Helac_Global.mod $(moddir)setscale.mod $(moddir)Helac_Func_1.mod $(moddir)Structf_PDFs.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)HO_chi2psia.o $(moddir)HO_chi2psia.mod : $(decaydir)HO_chi2psia.f90 $(moddir)Helac_ranmar.mod $(moddir)Kinetic_Func.mod $(moddir)HOVll.mod
	$(FC) -c $< $(FFLAGS)
makefile_pdf :
	(cd ./pdf ; make clean -f makefile_pdf ; make -f makefile_pdf ; cd ../)
makefile_fjcore :
	(cd ./jets/fjcore ; make clean -f makefile_fjcore ; make -f makefile_fjcore ; cd ../../)
makefile_fastjetfortran :
	(cd ./jets/fastjet ; make clean -f makefile_fastjetfortran ; make -f makefile_fastjetfortran ; cd ../../)
makefile_fjcore_heptoptagger :
	(cd ./analysis/heptoptagger ; make clean -f makefile_fjcore_heptoptagger ; make -f makefile_fjcore_heptoptagger ; cd ../../)
makefile_fastjet_heptoptagger :
	(cd ./analysis/heptoptagger ; make clean -f makefile_fastjet_heptoptagger ; make -f makefile_fastjet_heptoptagger ; cd ../../)
make_qedps :
	(cd ./shower/QEDPS; make clean ; make install; cp libqedps.a $(libdir)libqedps.a; rm *.a; make clean; cd ../../)
make_hbook :
	(cd ./analysis/hbook; make clean -f makefile_hbook ; make -f makefile_hbook; rm *.o; cd ../../)
make_minuit:
	(cd ./cernlib/minuit; make clean -f makefile_minuit ; make -f makefile_minuit; cd ../../)
clean:
	@(cd ./pdf ; make clean -f makefile_pdf ; cd ../)
	@(cd ./shower/QEDPS; make clean ; cd ../../)
	@(cd ./analysis/hbook; make clean -f makefile_hbook; cd ../../)
	@(cd ./cernlib/minuit; make clean -f makefile_minuit; cd ../../)
	@(cd ./analysis/HepMC; make clean_all -f makefile_hepmc2plot; cd ../../)
	@(rm -rf $(PWD)/addon/pp_psipsi_DPS/obj/*.o; rm -rf $(PWD)/addon/pp_psipsi_DPS/mod/*.mod; rm -rf $(PWD)/addon/pp_psipsi_DPS/bin/*; rm -rf $(PWD)/lib/*.a; rm -rf $(PWD)/addon/pp_psipsi_DPS/input/*~)
	@(rm -rf $(PWD)/addon/pp_psiX_CrystalBall/obj/*.o; rm -rf $(PWD)/addon/pp_psiX_CrystalBall/mod/*.mod; rm -rf $(PWD)/addon/pp_psiX_CrystalBall/bin/*; rm -rf $(PWD)/lib/*.a; rm -rf $(PWD)/addon/pp_psiX_CrystalBall/input/*~)
	@(rm -rf $(PWD)/addon/fit_pp_psiX_CrystalBall/obj/*.o; rm -rf $(PWD)/addon/fit_pp_psiX_CrystalBall/mod/*.mod; rm -rf $(PWD)/addon/fit_pp_psiX_CrystalBall/bin/*; rm -rf $(PWD)/lib/*.a; rm -rf $(PWD)/addon/fit_pp_psiX_CrystalBall/input/*~)
	@(rm -rf $(PWD)/addon/fit_pp_upsilonX_CrystalBall/obj/*.o; rm -rf $(PWD)/addon/fit_pp_upsilonX_CrystalBall/mod/*.mod; rm -rf $(PWD)/addon/fit_pp_upsilonX_CrystalBall/bin/*; rm -rf $(PWD)/lib/*.a; rm -rf $(PWD)/addon/fit_pp_upsilonX_CrystalBall/input/*~)
	@(rm -rf $(PWD)/addon/pp_aajj_DPS/bin/*; rm -rf $(PWD)/addon/pp_aajj_DPS/obj/*.o; rm -rf $(PWD)/addon/pp_aajj_DPS/mod/*.mod; rm -rf $(PWD)/addon/pp_aajj_DPS/input/*~)
	rm -rf $(obj) $(mod) $(exe) $(linkexe) $(PWD)/tmp/* *~ $(PWD)/input/*~ $(PWD)/src/*~ $(PWD)/output/* fort.* ./bin/ho_cluster ho_cluster ./bin/*
