HODIR=
FC=
objdir=$(HODIR)/addon/pp_psipsi_DPS/obj/
moddir=$(HODIR)/addon/pp_psipsi_DPS/mod/
srcdir=$(HODIR)/src/
thissrcdir=$(HODIR)/addon/pp_psipsi_DPS/src/
decaydir=$(HODIR)/decay/
libdir=$(HODIR)/lib/
showerdir=$(HODIR)/shower/interface/
lhapdfobj=
obj=$(objdir)Helac_Global.o $(objdir)Helac_Func_1.o  $(objdir)Constants.o $(objdir)CTEQ6PDF.o $(objdir)Kinetic_Func.o $(objdir)Structf_PDFs.o $(objdir)MC_VEGAS.o $(lhapdfobj) $(objdir)DecayInfo.o $(objdir)HOVll.o $(objdir)Decay_interface.o $(objdir)plot_pp_psipsi_dps.o $(objdir)pp_psipsi_DPS.o $(objdir)pp_psipsi_DPS_ME.o $(objdir)alfas_functions.o $(objdir)setscale.o $(objdir)Helac_ranmar.o $(objdir)QEDPS_interface.o $(objdir)pp_psipsi_dps_global.o $(objdir)Main_pp_psipsi_DPS.o $(objdir)Func_PSI.o $(objdir)HO_t2bw.o $(objdir)HO_chi2psia.o
mod = $(moddir)*.mod
exe = $(HODIR)/bin/HO_pp_psipsi_DPS
linkexe = $(HODIR)/addon/pp_psipsi_DPS/bin/HO_pp_psipsi_DPS
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
	$(FC) -o $(exe) $(obj) -L$(libdir) -lpdf -lgsdpdf -lhbook -lqedps -lminuit $(EXLIBS)

$(objdir)Main_pp_psipsi_DPS.o : $(thissrcdir)Main_pp_psipsi_DPS.f90 $(moddir)pp_psipsi_DPS.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)pp_psipsi_dps_global.o $(moddir)pp_psipsi_dps_global.mod: $(thissrcdir)pp_psipsi_dps_global.f90
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
$(objdir)DecayInfo.o $(moddir)DecayInfo.mod : $(decaydir)DecayInfo.f90 $(moddir)Helac_Global.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)HOVll.o $(moddir)HOVll.mod : $(decaydir)HOVll.f90 $(moddir)Helac_ranmar.mod $(moddir)Kinetic_Func.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)Decay_interface.o $(moddir)Decay_interface.mod : $(decaydir)Decay_interface.f90 $(moddir)Helac_Global.mod $(moddir)HOVll.mod $(moddir)Helac_Func_1.mod $(moddir)Kinetic_Func.mod $(moddir)QEDPS_interface.mod $(moddir)Constants.mod $(moddir)HO_t2bw.mod $(moddir)HO_chi2psia.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)QEDPS_interface.o $(moddir)QEDPS_interface.mod : $(showerdir)QEDPS_interface.f90 $(moddir)Helac_Global.mod $(moddir)Helac_Func_1.mod $(moddir)Kinetic_Func.mod $(moddir)setscale.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)pp_psipsi_DPS_ME.o $(moddir)pp_psipsi_DPS_ME.mod : $(thissrcdir)pp_psipsi_DPS_ME.f90 $(moddir)Helac_Global.mod $(moddir)pp_psipsi_dps_global.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)plot_pp_psipsi_dps.o $(moddir)plot_pp_psipsi_dps.mod : $(thissrcdir)plot_pp_psipsi_dps.f90 $(moddir)Helac_Global.mod $(moddir)Kinetic_Func.mod $(moddir)MC_VEGAS.mod $(moddir)Helac_Func_1.mod $(moddir)pp_psipsi_dps_global.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)pp_psipsi_DPS.o $(moddir)pp_psipsi_DPS.mod : $(thissrcdir)pp_psipsi_DPS.f90 $(moddir)Helac_Global.mod $(moddir)Helac_Func_1.mod $(moddir)Kinetic_Func.mod $(moddir)MC_VEGAS.mod $(moddir)Constants.mod $(moddir)Decay_interface.mod $(moddir)pp_psipsi_DPS_ME.mod $(moddir)plot_pp_psipsi_dps.mod $(moddir)DecayInfo.mod $(moddir)CTEQ6PDF.mod $(moddir)Structf_PDFs.mod $(moddir)HOVll.mod $(moddir)pp_psipsi_dps_global.mod $(moddir)Func_PSI.mod $(moddir)Helac_ranmar.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)pp_psipsi_DPS_LHAPDFs.o : $(thissrcdir)pp_psipsi_DPS_LHAPDFs.f90 $(moddir)Helac_Global.mod $(moddir)pp_psipsi_dps_global.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)Structf_LHAPDFs.o : $(srcdir)Structf_LHAPDFs.f90 $(moddir)Helac_Global.mod $(moddir)setscale.mod
	$(FC) -c $< $(FFLAGS)
$(objdir)HO_t2bw.o $(moddir)HO_t2bw.mod : $(decaydir)HO_t2bw.f90 $(moddir)Helac_ranmar.mod $(moddir)Kinetic_Func.mod
        $(FC) -c $< $(FFLAGS)
$(objdir)HO_chi2psia.o $(moddir)HO_chi2psia.mod : $(decaydir)HO_chi2psia.f90 $(moddir)Helac_ranmar.mod $(moddir)Kinetic_Func.mod $(moddir)HOVll.mod
	$(FC) -c $< $(FFLAGS)
makefile_gsdpdf:
	(cd ./pdf/gsdpdf ; make clean -f makefile_gsdpdf ; make -f makefile_gsdpdf ; cd ../)
clean:
	rm -rf $(obj) $(mod) $(exe) $(linkexe) $(HODIR)/tmp/* *~ $(HODIR)/input/*~ $(HODIR)/src/*~ $(thissrcdir)*~ $(HODIR)/output/* fort.*
end
