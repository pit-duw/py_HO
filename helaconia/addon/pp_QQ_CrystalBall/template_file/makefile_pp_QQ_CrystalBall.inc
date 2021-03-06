HODIR=
FC=
objdir=$(HODIR)/addon/pp_QQ_CrystalBall/obj/
moddir=$(HODIR)/addon/pp_QQ_CrystalBall/mod/
srcdir=$(HODIR)/src/
thissrcdir=$(HODIR)/addon/pp_QQ_CrystalBall/src/
decaydir=$(HODIR)/decay/
libdir=$(HODIR)/lib/
showerdir=$(HODIR)/shower/interface/
lhapdfobj=
obj=$(objdir)Helac_Global.o $(objdir)Helac_Func_1.o  $(objdir)Constants.o $(objdir)CTEQ6PDF.o $(objdir)Kinetic_Func.o $(objdir)Structf_PDFs.o $(objdir)MC_VEGAS.o $(lhapdfobj) $(objdir)plot_pp_QQ_cb.o $(objdir)pp_QQ_cb.o $(objdir)pp_QQ_cb_ME.o $(objdir)alfas_functions.o $(objdir)setscale.o $(objdir)Helac_ranmar.o $(objdir)pp_QQ_cb_global.o $(objdir)Main_pp_QQ_cb.o $(objdir)Func_PSI.o
mod = $(moddir)*.mod
exe = $(HODIR)/bin/HO_pp_QQ_CrystalBall
linkexe = $(HODIR)/addon/pp_QQ_CrystalBall/bin/HO_pp_QQ_CrystalBall
EXLIBS = 
all :   clean $(exe)
	ln -s $(exe) $(linkexe)

$(exe): $(obj)
	@cp *.o $(objdir)
	@cp *.mod $(moddir)
	@rm *.o *.mod
	$(FC) -o $(exe) $(obj) -L$(libdir) -lpdf -lhbook -lqedps -lminuit $(EXLIBS)

$(objdir)Main_pp_QQ_cb.o : $(thissrcdir)Main_pp_QQ_cb.f90 $(moddir)pp_QQ_cb.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)pp_QQ_cb_global.o $(moddir)pp_QQ_cb_global.mod: $(thissrcdir)pp_QQ_cb_global.f90
	$(FC) -c $< $(FFLAGS)
$(objdir)Helac_Global.o $(moddir)Helac_Global.mod : $(srcdir)Helac_Global.f90
	$(FC) -c $< $(FFLAGS)
$(objdir)Helac_Func_1.o $(moddir)Helac_Func_1.mod : $(srcdir)Helac_Func_1.f90 $(moddir)Helac_Global.mod $(moddir)alfas_functions.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)Constants.o $(moddir)Constants.mod: $(srcdir)Constants.f90 $(moddir)Helac_Global.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)Helac_ranmar.o $(moddir)Helac_ranmar.mod: $(srcdir)Helac_ranmar.f90
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
$(objdir)pp_QQ_cb_ME.o $(moddir)pp_QQ_cb_ME.mod : $(thissrcdir)pp_QQ_cb_ME.f90 $(moddir)Helac_Global.mod $(moddir)pp_QQ_cb_global.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)plot_pp_QQ_cb.o $(moddir)plot_pp_QQ_cb.mod : $(thissrcdir)plot_pp_QQ_cb.f90 $(moddir)Helac_Global.mod $(moddir)Kinetic_Func.mod $(moddir)MC_VEGAS.mod $(moddir)Helac_Func_1.mod $(moddir)pp_QQ_cb_global.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)pp_QQ_cb.o $(moddir)pp_QQ_cb.mod : $(thissrcdir)pp_QQ_cb.f90 $(moddir)Helac_Global.mod $(moddir)Helac_Func_1.mod $(moddir)Kinetic_Func.mod $(moddir)MC_VEGAS.mod $(moddir)Constants.mod $(moddir)pp_QQ_cb_ME.mod $(moddir)plot_pp_QQ_cb.mod $(moddir)CTEQ6PDF.mod $(moddir)Structf_PDFs.mod $(moddir)pp_QQ_cb_global.mod $(moddir)Func_PSI.mod $(moddir)Helac_ranmar.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)pp_QQ_cb_LHAPDFs.o : $(thissrcdir)pp_QQ_cb_LHAPDFs.f90 $(moddir)Helac_Global.mod $(moddir)pp_QQ_cb_global.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)Structf_LHAPDFs.o : $(srcdir)Structf_LHAPDFs.f90 $(moddir)Helac_Global.mod $(moddir)setscale.mod
	$(FC) -c $< $(FFLAGS)
clean:
	rm -rf $(obj) $(mod) $(exe) $(linkexe) $(HODIR)/tmp/* *~ $(HODIR)/input/*~ $(HODIR)/src/*~ $(thissrcdir)*~ $(HODIR)/output/* fort.*
end
