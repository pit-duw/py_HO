HODIR=
FC=
objdir=$(HODIR)/addon/fit_pp_psiX_CrystalBall/obj/
moddir=$(HODIR)/addon/fit_pp_psiX_CrystalBall/mod/
srcdir=$(HODIR)/src/
thissrcdir=$(HODIR)/addon/fit_pp_psiX_CrystalBall/src/
libdir=$(HODIR)/lib/
lhapdfobj=
obj=$(objdir)plot_fit_cb.o $(objdir)fit_cb_global.o $(objdir)fit_crystalball.o $(objdir)Constants.o $(objdir)Helac_Global.o $(objdir)MC_VEGAS.o $(objdir)Func_PSI.o $(objdir)CTEQ6PDF.o $(objdir)Structf_PDFs.o $(objdir)Kinetic_Func.o $(objdir)Helac_Func_1.o $(objdir)setscale.o $(objdir)alfas_functions.o $(lhapdfobj)
mod = $(moddir)*.mod
exe = $(HODIR)/bin/HO_fit_pp_psiX_CrystalBall
linkexe = $(HODIR)/addon/fit_pp_psiX_CrystalBall/bin/HO_fit_pp_psiX_CrystalBall
EXLIBS =
all :   clean $(exe)
	ln -s $(exe) $(linkexe)

$(exe): $(obj)
	@cp *.o $(objdir)
	@cp *.mod $(moddir)
	@rm *.o *.mod
	$(FC) -o $(exe) $(obj) -L$(libdir) -lpdf -lminuit $(EXLIBS)

$(objdir)fit_cb_global.o $(moddir)fit_cb_global.mod: $(thissrcdir)fit_cb_global.f90
	$(FC) -c $< $(FFLAGS)
$(objdir)fit_crystalball.o: $(thissrcdir)fit_crystalball.f90 $(moddir)fit_cb_global.mod $(moddir)Constants.mod $(moddir)Helac_Global.mod $(moddir)MC_VEGAS.mod $(moddir)Func_PSI.mod $(moddir)CTEQ6PDF.mod $(moddir)Structf_PDFs.mod $(moddir)Kinetic_Func.mod $(moddir)plot_fit_cb.mod
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
$(objdir)setscale.o $(moddir)setscale.mod: $(srcdir)setscale.f90
	$(FC) -c $< $(FFLAGS)
$(objdir)CTEQ6PDF.o $(moddir)CTEQ6PDF.mod: $(srcdir)CTEQ6PDF.f90
	$(FC) -c $< $(FFLAGS)
$(objdir)Kinetic_Func.o $(moddir)Kinetic_Func.mod: $(srcdir)Kinetic_Func.f90 $(moddir)Helac_Global.mod $(moddir)Helac_Func_1.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)Structf_LHAPDFs.o : $(srcdir)Structf_LHAPDFs.f90 $(moddir)Helac_Global.mod $(moddir)setscale.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)MC_VEGAS.o $(moddir)MC_VEGAS.mod: $(srcdir)MC_VEGAS.f90
	$(FC) -c $< $(FFLAGS)
$(objdir)alfas_functions.o $(moddir)alfas_functions.mod: $(srcdir)alfas_functions.f90 $(moddir)Helac_Global.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)fit_cb_LHAPDFs.o : $(thissrcdir)fit_cb_LHAPDFs.f90 $(moddir)fit_cb_global.mod $(moddir)Helac_Global.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)plot_fit_cb.o $(moddir)plot_fit_cb.mod: $(thissrcdir)plot_fit_cb.f90
	$(FC) -c $< $(FFLAGS)
clean:
	rm -rf $(obj) $(mod) $(exe) $(linkexe) $(HODIR)/tmp/* *~ $(HODIR)/input/*~ $(HODIR)/src/*~ $(thissrcdir)*~ $(HODIR)/output/* fort.*
end
