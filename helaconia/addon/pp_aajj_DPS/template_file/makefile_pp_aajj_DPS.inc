HODIR=
FC=
objdir=$(HODIR)/addon/pp_aajj_DPS/obj/
moddir=$(HODIR)/addon/pp_aajj_DPS/mod/
srcdir=$(HODIR)/src/
thissrcdir=$(HODIR)/addon/pp_aajj_DPS/src/
libdir=$(HODIR)/lib/
showerdir=$(HODIR)/shower/interface/
lhapdfobj=
obj=$(objdir)Helac_Global.o $(objdir)Helac_Func_1.o  $(objdir)Constants.o $(objdir)CTEQ6PDF.o $(objdir)Kinetic_Func.o $(objdir)Structf_PDFs.o $(objdir)MC_VEGAS.o $(lhapdfobj) $(objdir)pp_aajj_dps_func.o $(objdir)plot_pp_aajj_dps.o $(objdir)pp_aajj_DPS.o $(objdir)Amp_gg_aa.o $(objdir)ME2_pp_jj.o $(objdir)ME2_pp_aj.o $(objdir)ME2_pp_aa.o $(objdir)alfas_functions.o $(objdir)setscale.o $(objdir)pp_aajj_dps_global.o $(objdir)Main_pp_aajj_DPS.o $(objdir)Func_PSI.o $(objdir)pp_aajj_DPS_cuts.o
mod = $(moddir)*.mod
exe = $(HODIR)/bin/HO_pp_aajj_DPS
linkexe = $(HODIR)/addon/pp_aajj_DPS/bin/HO_pp_aajj_DPS
EXLIBS = 

ifneq ("$(wildcard $(HODIR)/lib/libgsdpdf.a)","")
      makegsdpdf = 
else
      makegsdpdf = makefile_gsdpdf
endif

all :   clean $(exe)
	ln -s $(exe) $(linkexe)

$(exe): $(makegsdpdf) $(obj)
	@cp *.o $(objdir)
	@cp *.mod $(moddir)
	@rm *.o *.mod
	$(FC) -o $(exe) $(obj) -L$(libdir) -lpdf -lgsdpdf -lhbook $(EXLIBS)

$(objdir)Main_pp_aajj_DPS.o : $(thissrcdir)Main_pp_aajj_DPS.f90 $(moddir)pp_aajj_DPS.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)pp_aajj_dps_global.o $(moddir)pp_aajj_dps_global.mod: $(thissrcdir)pp_aajj_dps_global.f90
	$(FC) -c $< $(FFLAGS)
$(objdir)Helac_Global.o $(moddir)Helac_Global.mod : $(srcdir)Helac_Global.f90
	$(FC) -c $< $(FFLAGS)
$(objdir)Helac_Func_1.o $(moddir)Helac_Func_1.mod : $(srcdir)Helac_Func_1.f90 $(moddir)Helac_Global.mod $(moddir)alfas_functions.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)Constants.o $(moddir)Constants.mod: $(srcdir)Constants.f90 $(moddir)Helac_Global.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)Func_PSI.o $(moddir)Func_PSI.mod: $(srcdir)Func_PSI.f90
	$(FC) -c $< $(FFLAGS)			   
$(objdir)Structf_PDFs.o $(moddir)Structf_PDFs.mod:$(srcdir)Structf_PDFs.f90 $(moddir)Helac_Global.mod $(moddir)CTEQ6PDF.mod $(moddir)setscale.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)CTEQ6PDF.o $(moddir)CTEQ6PDF.mod: $(srcdir)CTEQ6PDF.f90
	$(FC) -c $< $(FFLAGS)
$(objdir)Kinetic_Func.o $(moddir)Kinetic_Func.mod: $(srcdir)Kinetic_Func.f90 $(moddir)Helac_Global.mod $(moddir)Helac_Func_1.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)MC_VEGAS.o $(moddir)MC_VEGAS.mod: $(srcdir)MC_VEGAS.f90
	$(FC) -c $< $(FFLAGS)
$(objdir)alfas_functions.o $(moddir)alfas_functions.mod: $(srcdir)alfas_functions.f90 $(moddir)Helac_Global.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)setscale.o $(moddir)setscale.mod: $(srcdir)setscale.f90
	$(FC) -c $< $(FFLAGS)
$(objdir)ME2_pp_aa.o $(moddir)ME2_pp_aa.mod : $(thissrcdir)ME2_pp_aa.f90 $(moddir)Helac_Global.mod $(moddir)pp_aajj_dps_global.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)ME2_pp_aj.o $(moddir)ME2_pp_aj.mod : $(thissrcdir)ME2_pp_aj.f90 $(moddir)Helac_Global.mod $(moddir)pp_aajj_dps_global.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)ME2_pp_jj.o $(moddir)ME2_pp_jj.mod : $(thissrcdir)ME2_pp_jj.f90 $(moddir)Helac_Global.mod $(moddir)pp_aajj_dps_global.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)Amp_gg_aa.o : $(thissrcdir)Amp_gg_aa.f
	$(FC) -c $< $(FFLAGS)
$(objdir)pp_aajj_dps_func.o $(moddir)pp_aajj_dps_func.mod : $(thissrcdir)pp_aajj_dps_func.f90 $(moddir)pp_aajj_dps_global.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)plot_pp_aajj_dps.o $(moddir)plot_pp_aajj_dps.mod : $(thissrcdir)plot_pp_aajj_dps.f90 $(moddir)Helac_Global.mod $(moddir)Kinetic_Func.mod $(moddir)MC_VEGAS.mod $(moddir)Helac_Func_1.mod $(moddir)pp_aajj_dps_global.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)pp_aajj_DPS.o $(moddir)pp_aajj_DPS.mod : $(thissrcdir)pp_aajj_DPS.f90 $(moddir)Helac_Global.mod $(moddir)Helac_Func_1.mod $(moddir)Kinetic_Func.mod $(moddir)MC_VEGAS.mod $(moddir)Constants.mod $(moddir)plot_pp_aajj_dps.mod $(moddir)CTEQ6PDF.mod $(moddir)Structf_PDFs.mod $(moddir)pp_aajj_dps_global.mod $(moddir)Func_PSI.mod $(moddir)ME2_pp_jj.mod $(moddir)ME2_pp_aj.mod $(moddir)ME2_pp_aa.mod $(moddir)pp_aajj_dps_func.mod $(moddir)pp_aajj_DPS_cuts.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)pp_aajj_DPS_LHAPDFs.o : $(thissrcdir)pp_aajj_DPS_LHAPDFs.f90 $(moddir)Helac_Global.mod $(moddir)pp_aajj_dps_global.mod $(moddir)pp_aajj_dps_func.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)Structf_LHAPDFs.o : $(srcdir)Structf_LHAPDFs.f90 $(moddir)Helac_Global.mod $(moddir)setscale.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)pp_aajj_DPS_cuts.o $(moddir)pp_aajj_DPS_cuts.mod : $(thissrcdir)pp_aajj_DPS_cuts.f90 $(moddir)Helac_Global.mod $(moddir)pp_aajj_dps_global.mod $(moddir)Constants.mod $(moddir)Kinetic_Func.mod $(moddir)MC_VEGAS.mod
	$(FC) -c $< $(FFLAGS)
makefile_gsdpdf:
	(cd ./pdf/gsdpdf ; make clean -f makefile_gsdpdf ; make -f makefile_gsdpdf ; cd ../)
clean:
	rm -rf $(obj) $(mod) $(exe) $(linkexe) $(HODIR)/tmp/* *~ $(HODIR)/input/*~ $(HODIR)/src/*~ $(thissrcdir)*~ $(HODIR)/output/* fort.*
end
