
  reset
  set lmargin 10
  set rmargin 0
  set terminal postscript col enhanced
  set output "fit_results.ps"
  set style data histeps
  set key off
  set style line 1 lt 1 lc rgb "#006D4F" lw 1.8
  set style line 2 lt 1 lc rgb "#B90091" lw 1.8

  set label "{/Symbol @\326\140}S=   7.00TeV,  2.00<y<  2.50" font ",14" at graph 0.1, graph 0.94
  set title "SET NUMBER =  1 distribution" font "Helvetica, 20"
  set xlabel "P_T [GeV]" font "Helvetica, 20"
  set ylabel "d^2{/Symbol s}/dP_Tdy [nb/GeV]" font "Helvetica, 20"
  set label front "HELAC-ONIA" font "Helvetica, 15" rotate by 90 at graph 1.02, graph 0.4
  set xrange [    0.00000:  14.00000]
  set format y "10^{%T}"
  set key horizontal noreverse maxcols 1 width -4 at graph 0.92, graph 0.9
  set logscale y
  plot \
  "./comparison_1.dat" using 1:($4):($4+$5):($4-$5) with yerror ls 1 title "EX",\
  "./comparison_1.dat" using 1:($2) with histeps ls 2 title "TH"
  unset label

  set label "{/Symbol @\326\140}S=   7.00TeV,  2.50<y<  3.00" font ",14" at graph 0.1, graph 0.94
  set title "SET NUMBER =  2 distribution" font "Helvetica, 20"
  set xlabel "P_T [GeV]" font "Helvetica, 20"
  set ylabel "d^2{/Symbol s}/dP_Tdy [nb/GeV]" font "Helvetica, 20"
  set label front "HELAC-ONIA" font "Helvetica, 15" rotate by 90 at graph 1.02, graph 0.4
  set xrange [    0.00000:  14.00000]
  set format y "10^{%T}"
  set key horizontal noreverse maxcols 1 width -4 at graph 0.92, graph 0.9
  set logscale y
  plot \
  "./comparison_2.dat" using 1:($4):($4+$5):($4-$5) with yerror ls 1 title "EX",\
  "./comparison_2.dat" using 1:($2) with histeps ls 2 title "TH"
  unset label

  set label "{/Symbol @\326\140}S=   7.00TeV,  3.00<y<  3.50" font ",14" at graph 0.1, graph 0.94
  set title "SET NUMBER =  3 distribution" font "Helvetica, 20"
  set xlabel "P_T [GeV]" font "Helvetica, 20"
  set ylabel "d^2{/Symbol s}/dP_Tdy [nb/GeV]" font "Helvetica, 20"
  set label front "HELAC-ONIA" font "Helvetica, 15" rotate by 90 at graph 1.02, graph 0.4
  set xrange [    0.00000:  14.00000]
  set format y "10^{%T}"
  set key horizontal noreverse maxcols 1 width -4 at graph 0.92, graph 0.9
  set logscale y
  plot \
  "./comparison_3.dat" using 1:($4):($4+$5):($4-$5) with yerror ls 1 title "EX",\
  "./comparison_3.dat" using 1:($2) with histeps ls 2 title "TH"
  unset label

  set label "{/Symbol @\326\140}S=   7.00TeV,  3.50<y<  4.00" font ",14" at graph 0.1, graph 0.94
  set title "SET NUMBER =  4 distribution" font "Helvetica, 20"
  set xlabel "P_T [GeV]" font "Helvetica, 20"
  set ylabel "d^2{/Symbol s}/dP_Tdy [nb/GeV]" font "Helvetica, 20"
  set label front "HELAC-ONIA" font "Helvetica, 15" rotate by 90 at graph 1.02, graph 0.4
  set xrange [    0.00000:  13.00000]
  set format y "10^{%T}"
  set key horizontal noreverse maxcols 1 width -4 at graph 0.92, graph 0.9
  set logscale y
  plot \
  "./comparison_4.dat" using 1:($4):($4+$5):($4-$5) with yerror ls 1 title "EX",\
  "./comparison_4.dat" using 1:($2) with histeps ls 2 title "TH"
  unset label

  set label "{/Symbol @\326\140}S=   7.00TeV,  4.00<y<  4.50" font ",14" at graph 0.1, graph 0.94
  set title "SET NUMBER =  5 distribution" font "Helvetica, 20"
  set xlabel "P_T [GeV]" font "Helvetica, 20"
  set ylabel "d^2{/Symbol s}/dP_Tdy [nb/GeV]" font "Helvetica, 20"
  set label front "HELAC-ONIA" font "Helvetica, 15" rotate by 90 at graph 1.02, graph 0.4
  set xrange [    0.00000:  11.00000]
  set format y "10^{%T}"
  set key horizontal noreverse maxcols 1 width -4 at graph 0.92, graph 0.9
  set logscale y
  plot \
  "./comparison_5.dat" using 1:($4):($4+$5):($4-$5) with yerror ls 1 title "EX",\
  "./comparison_5.dat" using 1:($2) with histeps ls 2 title "TH"
  unset label

  set label "{/Symbol @\326\140}S=   7.00TeV,  0.00<y<  0.75" font ",14" at graph 0.1, graph 0.94
  set title "SET NUMBER =  6 distribution" font "Helvetica, 20"
  set xlabel "P_T [GeV]" font "Helvetica, 20"
  set ylabel "d^2{/Symbol s}/dP_Tdy [nb/GeV]" font "Helvetica, 20"
  set label front "HELAC-ONIA" font "Helvetica, 15" rotate by 90 at graph 1.02, graph 0.4
  set xrange [    7.00000:  70.00000]
  set format y "10^{%T}"
  set key horizontal noreverse maxcols 1 width -4 at graph 0.92, graph 0.9
  set logscale y
  plot \
  "./comparison_6.dat" using 1:($4):($4+$5):($4-$5) with yerror ls 1 title "EX",\
  "./comparison_6.dat" using 1:($2) with histeps ls 2 title "TH"
  unset label

  set label "{/Symbol @\326\140}S=   7.00TeV,  0.75<y<  1.50" font ",14" at graph 0.1, graph 0.94
  set title "SET NUMBER =  7 distribution" font "Helvetica, 20"
  set xlabel "P_T [GeV]" font "Helvetica, 20"
  set ylabel "d^2{/Symbol s}/dP_Tdy [nb/GeV]" font "Helvetica, 20"
  set label front "HELAC-ONIA" font "Helvetica, 15" rotate by 90 at graph 1.02, graph 0.4
  set xrange [    5.00000:  70.00000]
  set format y "10^{%T}"
  set key horizontal noreverse maxcols 1 width -4 at graph 0.92, graph 0.9
  set logscale y
  plot \
  "./comparison_7.dat" using 1:($4):($4+$5):($4-$5) with yerror ls 1 title "EX",\
  "./comparison_7.dat" using 1:($2) with histeps ls 2 title "TH"
  unset label

  set label "{/Symbol @\326\140}S=   7.00TeV,  1.50<y<  2.00" font ",14" at graph 0.1, graph 0.94
  set title "SET NUMBER =  8 distribution" font "Helvetica, 20"
  set xlabel "P_T [GeV]" font "Helvetica, 20"
  set ylabel "d^2{/Symbol s}/dP_Tdy [nb/GeV]" font "Helvetica, 20"
  set label front "HELAC-ONIA" font "Helvetica, 15" rotate by 90 at graph 1.02, graph 0.4
  set xrange [    1.00000:  30.00000]
  set format y "10^{%T}"
  set key horizontal noreverse maxcols 1 width -4 at graph 0.92, graph 0.9
  set logscale y
  plot \
  "./comparison_8.dat" using 1:($4):($4+$5):($4-$5) with yerror ls 1 title "EX",\
  "./comparison_8.dat" using 1:($2) with histeps ls 2 title "TH"
  unset label

  set label "{/Symbol @\326\140}S=   7.00TeV,  2.00<y<  2.40" font ",14" at graph 0.1, graph 0.94
  set title "SET NUMBER =  9 distribution" font "Helvetica, 20"
  set xlabel "P_T [GeV]" font "Helvetica, 20"
  set ylabel "d^2{/Symbol s}/dP_Tdy [nb/GeV]" font "Helvetica, 20"
  set label front "HELAC-ONIA" font "Helvetica, 15" rotate by 90 at graph 1.02, graph 0.4
  set xrange [    5.00000:  30.00000]
  set format y "10^{%T}"
  set key horizontal noreverse maxcols 1 width -4 at graph 0.92, graph 0.9
  set logscale y
  plot \
  "./comparison_9.dat" using 1:($4):($4+$5):($4-$5) with yerror ls 1 title "EX",\
  "./comparison_9.dat" using 1:($2) with histeps ls 2 title "TH"
  unset label

  set label "{/Symbol @\326\140}S=   7.00TeV,  0.00<y<  0.90" font ",14" at graph 0.1, graph 0.94
  set title "SET NUMBER = 10 distribution" font "Helvetica, 20"
  set xlabel "P_T [GeV]" font "Helvetica, 20"
  set ylabel "d^2{/Symbol s}/dP_Tdy [nb/GeV]" font "Helvetica, 20"
  set label front "HELAC-ONIA" font "Helvetica, 15" rotate by 90 at graph 1.02, graph 0.4
  set xrange [    8.00000:  70.00000]
  set format y "10^{%T}"
  set key horizontal noreverse maxcols 1 width -4 at graph 0.92, graph 0.9
  set logscale y
  plot \
  "./comparison_10.dat" using 1:($4):($4+$5):($4-$5) with yerror ls 1 title "EX",\
  "./comparison_10.dat" using 1:($2) with histeps ls 2 title "TH"
  unset label

  set label "{/Symbol @\326\140}S=   7.00TeV,  0.90<y<  1.20" font ",14" at graph 0.1, graph 0.94
  set title "SET NUMBER = 11 distribution" font "Helvetica, 20"
  set xlabel "P_T [GeV]" font "Helvetica, 20"
  set ylabel "d^2{/Symbol s}/dP_Tdy [nb/GeV]" font "Helvetica, 20"
  set label front "HELAC-ONIA" font "Helvetica, 15" rotate by 90 at graph 1.02, graph 0.4
  set xrange [    8.00000:  45.00000]
  set format y "10^{%T}"
  set key horizontal noreverse maxcols 1 width -4 at graph 0.92, graph 0.9
  set logscale y
  plot \
  "./comparison_11.dat" using 1:($4):($4+$5):($4-$5) with yerror ls 1 title "EX",\
  "./comparison_11.dat" using 1:($2) with histeps ls 2 title "TH"
  unset label

  set label "{/Symbol @\326\140}S=   7.00TeV,  1.20<y<  1.60" font ",14" at graph 0.1, graph 0.94
  set title "SET NUMBER = 12 distribution" font "Helvetica, 20"
  set xlabel "P_T [GeV]" font "Helvetica, 20"
  set ylabel "d^2{/Symbol s}/dP_Tdy [nb/GeV]" font "Helvetica, 20"
  set label front "HELAC-ONIA" font "Helvetica, 15" rotate by 90 at graph 1.02, graph 0.4
  set xrange [    6.50000:  45.00000]
  set format y "10^{%T}"
  set key horizontal noreverse maxcols 1 width -4 at graph 0.92, graph 0.9
  set logscale y
  plot \
  "./comparison_12.dat" using 1:($4):($4+$5):($4-$5) with yerror ls 1 title "EX",\
  "./comparison_12.dat" using 1:($2) with histeps ls 2 title "TH"
  unset label

  set label "{/Symbol @\326\140}S=   7.00TeV,  1.60<y<  2.10" font ",14" at graph 0.1, graph 0.94
  set title "SET NUMBER = 13 distribution" font "Helvetica, 20"
  set xlabel "P_T [GeV]" font "Helvetica, 20"
  set ylabel "d^2{/Symbol s}/dP_Tdy [nb/GeV]" font "Helvetica, 20"
  set label front "HELAC-ONIA" font "Helvetica, 15" rotate by 90 at graph 1.02, graph 0.4
  set xrange [    6.50000:  30.00000]
  set format y "10^{%T}"
  set key horizontal noreverse maxcols 1 width -4 at graph 0.92, graph 0.9
  set logscale y
  plot \
  "./comparison_13.dat" using 1:($4):($4+$5):($4-$5) with yerror ls 1 title "EX",\
  "./comparison_13.dat" using 1:($2) with histeps ls 2 title "TH"
  unset label

  set label "{/Symbol @\326\140}S=   7.00TeV,  2.10<y<  2.40" font ",14" at graph 0.1, graph 0.94
  set title "SET NUMBER = 14 distribution" font "Helvetica, 20"
  set xlabel "P_T [GeV]" font "Helvetica, 20"
  set ylabel "d^2{/Symbol s}/dP_Tdy [nb/GeV]" font "Helvetica, 20"
  set label front "HELAC-ONIA" font "Helvetica, 15" rotate by 90 at graph 1.02, graph 0.4
  set xrange [    5.50000:  30.00000]
  set format y "10^{%T}"
  set key horizontal noreverse maxcols 1 width -4 at graph 0.92, graph 0.9
  set logscale y
  plot \
  "./comparison_14.dat" using 1:($4):($4+$5):($4-$5) with yerror ls 1 title "EX",\
  "./comparison_14.dat" using 1:($2) with histeps ls 2 title "TH"
